#!/usr/bin/env python3
import argparse
import importlib
import json
import multiprocessing as mp
import os
import re
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent

MODEL_IDS = {
    "HG chipseq": "1VAsp8I904y_J0PUhAQqpSlCn1IqfG0FB",
    "HG kouzine": "1dAeAt5Gu2cadwDhbc7OnenUgDLHlUvkx",
    "MM curax":   "1W6GEgHNoitlB-xXJbLJ_jDW4BF35W1Sd",
    "MM kouzine": "1dXpQFmheClKXIEoqcZ7kgCwx6hzVCv3H",
}
AUX_IDS = {
    "config.json":             "10sF8Ywktd96HqAL0CwvlZZUUGj05CGk5",
    "special_tokens_map.json": "16bT7HDv71aRwyh3gBUbKwign1mtyLD2d",
    "tokenizer_config.json":   "1EE9goZ2JRSD8UTx501q71lGCk-CK3kqG",
    "vocab.txt":               "1gZZdtAoDnDiLQqjQfGyuwt268Pe5sXW0",
}


def download_model(model_dir: Path, model_name: str):
    """Скачать веса и файлы токенизатора с Google Drive (только недостающие)."""
    model_dir.mkdir(parents=True, exist_ok=True)
    need = {}
    if not (model_dir / "pytorch_model.bin").exists():
        need["pytorch_model.bin"] = MODEL_IDS[model_name]
    for fname, fid in AUX_IDS.items():
        if not (model_dir / fname).exists():
            need[fname] = fid
    if not need:
        print("модель уже на месте:", sorted(p.name for p in model_dir.iterdir()))
        return
    import gdown  # импортируем только если реально надо качать
    for fname, fid in need.items():
        gdown.download(id=fid, output=str(model_dir / fname), quiet=False)
    print("файлы модели:", sorted(p.name for p in model_dir.iterdir()))


def merge_intervals(intervals):
    """Слить пересекающиеся/смежные (start, end) — убирает дубли из перекрытия кусков."""
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [list(intervals[0])]
    for s, e in intervals[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return merged


def _safe(name):
    return re.sub(r"[^\w.-]", "_", str(name))


def build_tasks(fasta_in, chroms, max_bp, chunk, overlap, threshold, min_len,
                batch_size, cache_dir):
    """Нарезать геном на куски-задачи. id куска = {chrom}_{start}_{end}."""
    from Bio import SeqIO
    tasks = []
    for rec in SeqIO.parse(str(fasta_in), "fasta"):
        if chroms is not None and rec.id not in chroms:
            continue
        seq = str(rec.seq)
        L = len(seq) if max_bp is None else min(len(seq), max_bp)
        for st in range(0, L, chunk):
            en = min(st + chunk + overlap, L)
            cache_path = str(cache_dir / f"{_safe(rec.id)}_{st}_{en}.tsv")
            tasks.append((rec.id, st, seq[st:en], threshold, min_len, batch_size, cache_path))
    return tasks


def check_or_write_params(cache_dir, model_dir, fasta_in, chunk, overlap, threshold, min_len):
    """Сверить параметры, влияющие на результат, с сохранёнными в кэше (или записать)."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    st = Path(fasta_in).stat()
    current = {
        "model_dir": str(model_dir),
        "genome": str(fasta_in), "genome_size": st.st_size, "genome_mtime": int(st.st_mtime),
        "CHUNK": chunk, "OVERLAP": overlap,
        "THRESHOLD": threshold, "MIN_LEN": min_len,
    }
    pfile = cache_dir / "params.json"
    if pfile.exists():
        saved = json.loads(pfile.read_text())
        diff = {k: (saved.get(k), v) for k, v in current.items() if saved.get(k) != v}
        if diff:
            raise ValueError(f"Параметры кэша {cache_dir} не совпадают с текущими: {diff}. "
                             "Очисти папку кэша или укажи другую CACHE_DIR.")
        print("кэш: параметры совпали ->", pfile)
    else:
        pfile.write_text(json.dumps(current, indent=2))
        print("кэш: параметры сохранены ->", pfile)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fasta", type=Path, default=BASE_DIR / "data" / "genome.fa")
    p.add_argument("--bed-out", type=Path, default=BASE_DIR / "results" / "bed" / "zdnabert.bed")
    p.add_argument("--cache-dir", type=Path, default=BASE_DIR / "results" / "cache" / "zdnabert")
    p.add_argument("--model-dir", type=Path, default=BASE_DIR / "models" / "6-new-12w-0")
    p.add_argument("--model", default="HG kouzine", choices=list(MODEL_IDS))
    p.add_argument("--chroms", nargs="*", default=["1"],
                   help="имена хромосом; пусто или 'all' = все")
    p.add_argument("--max-bp", default="2400000",
                   help="bp на хромосому; 'None'/'all' = вся хромосома")
    p.add_argument("--chunk", type=int, default=200_000)
    p.add_argument("--overlap", type=int, default=1_000)
    p.add_argument("--threshold", type=float, default=0.5)
    p.add_argument("--min-len", type=int, default=10)
    p.add_argument("--n-proc", type=int, default=6)
    p.add_argument("--threads-per", type=int, default=16)
    p.add_argument("--batch-size", type=int, default=64)
    p.add_argument("--skip-download", action="store_true",
                   help="не скачивать модель (считать, что уже есть)")
    return p.parse_args()


def main():
    args = parse_args()

    chroms = None if (not args.chroms or args.chroms == ["all"]) else args.chroms
    max_bp = None if str(args.max_bp).lower() in ("none", "all") else int(args.max_bp)

    args.bed_out.parent.mkdir(parents=True, exist_ok=True)

    if not args.skip_download:
        download_model(args.model_dir, args.model)
    if not (args.model_dir / "pytorch_model.bin").exists():
        raise FileNotFoundError(f"нет модели в {args.model_dir} — убери --skip-download")
    if not args.fasta.exists():
        raise FileNotFoundError(f"нет генома {args.fasta} — сначала скачай и распакуй геном")

    # чтобы spawn-воркеры нашли модуль zdna_worker
    sys.path.insert(0, str(BASE_DIR))
    import zdna_worker as zw
    importlib.reload(zw)

    check_or_write_params(args.cache_dir, args.model_dir, args.fasta,
                          args.chunk, args.overlap, args.threshold, args.min_len)

    tasks = build_tasks(args.fasta, chroms, max_bp, args.chunk, args.overlap,
                        args.threshold, args.min_len, args.batch_size, args.cache_dir)
    cached = [t for t in tasks if os.path.exists(t[-1])]
    todo = [t for t in tasks if not os.path.exists(t[-1])]
    print(f"кусков: {len(tasks)} | из кэша: {len(cached)} | считать: {len(todo)} "
          f"| параллель: {args.n_proc} x {args.threads_per}")

    try:
        from tqdm import tqdm
    except ImportError:
        def tqdm(it, **kw):
            return it

    results = []
    # кэш-хиты читаем здесь — без пула и без загрузки модели
    for t in tqdm(cached, total=len(cached), desc="кэш"):
        results.extend(zw.process_chunk(t))

    # остальное считаем параллельно (пул только если есть что считать)
    if todo:
        ctx = mp.get_context("spawn")
        with ProcessPoolExecutor(max_workers=args.n_proc, mp_context=ctx,
                                 initializer=zw.init_worker,
                                 initargs=(str(args.model_dir), args.threads_per)) as ex:
            futures = [ex.submit(zw.process_chunk, t) for t in todo]
            for f in tqdm(as_completed(futures), total=len(futures), desc="куски"):
                results.extend(f.result())
    else:
        print("всё из кэша — инференс не нужен")

    # сливаем интервалы по хромосомам и пишем BED
    by_chrom = defaultdict(list)
    for chrom, s, e in results:
        by_chrom[chrom].append((s, e))

    n_written = 0
    with open(args.bed_out, "w") as bed:
        for chrom in sorted(by_chrom):
            for s, e in merge_intervals(by_chrom[chrom]):
                bed.write(f"{chrom}\t{s}\t{e}\tZDNABERT\t{e - s}\t.\n")
                n_written += 1
    print("готово ->", args.bed_out, "| интервалов:", n_written)


if __name__ == "__main__":
    main()
