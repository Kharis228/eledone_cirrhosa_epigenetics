#!/usr/bin/env python3
"""Вторичные структуры ДНК у Eledone cirrhosa (пункт 4) — скрипт-версия ноутбука.

Делает два BED по всему геному (или подмножеству хромосом):
  - results/bed/g4.bed    — G-квадруплексы (regex на обоих стрендах, верхний регистр)
  - results/bed/zhunt.bed — Z-ДНК (zhunt3, z-score > порога)

zhunt однопоточный, поэтому каждая хромосома режется на куски и считается
параллельно (ThreadPool + внешний zhunt3 + awk-фильтр). Геном обрабатывается
потоково по одной хромосоме за раз, чтобы не держать все 3 Гб в памяти.

Зависимости: biopython; внешние gcc, wget, awk (стандартно в Linux).

Примеры:
    python structures_predict.py                       # весь геном
    python structures_predict.py --chroms 1 2          # только chr1, chr2
    python structures_predict.py --only g4             # только квадруплексы
    python structures_predict.py --zhunt-proc 48 --z-threshold 500
"""

import argparse
import os
import re
import subprocess
import tempfile
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

try:
    from tqdm import tqdm
except ImportError:                       # без tqdm — просто без прогресс-бара
    def tqdm(it, **kw):
        return it

BASE_DIR = Path(__file__).resolve().parent

# G-квадруплексы (оба стренда, верхний регистр)
G4_PLUS = re.compile(r"(?:G{3,5}[ACGT]{1,7}){3,}G{3,5}")
G4_MINUS = re.compile(r"(?:C{3,5}[ACGT]{1,7}){3,}C{3,5}")

ZHUNT_SRC_URL = (
    "https://raw.githubusercontent.com/diegopenilla/FlaskHunt/master/zhunt3.c"
)


def build_zhunt(scripts_dir: Path) -> Path:
    """Скачать (если надо) и скомпилировать zhunt3. Вернуть путь к бинарю."""
    scripts_dir.mkdir(parents=True, exist_ok=True)
    src = scripts_dir / "zhunt3.c"
    binary = scripts_dir / "zhunt3"
    if not src.exists():
        subprocess.run(["wget", "-q", "-O", str(src), ZHUNT_SRC_URL], check=True)
    subprocess.run(["gcc", str(src), "-o", str(binary), "-lm"], check=True)
    return binary


def list_chrom_ids(fasta, chroms):
    """Список id хромосом, которые будут обработаны (читает только заголовки FASTA)."""
    ids = []
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                rid = line[1:].split()[0]
                if chroms is None or rid in chroms:
                    ids.append(rid)
    return ids


def find_g4(seq, chrom):
    rows = []
    for strand, pat in (("+", G4_PLUS), ("-", G4_MINUS)):
        for m in pat.finditer(seq):
            rows.append((chrom, m.start(), m.end(), strand))
    rows.sort(key=lambda r: (r[1], r[2]))
    return rows


def merge_positions(positions):
    """Слить последовательные геномные позиции в интервалы (start, end)."""
    positions = sorted(positions)
    out = []
    if positions:
        st = pr = positions[0]
        for v in positions[1:]:
            if v <= pr + 1:
                pr = v
            else:
                out.append((st, pr + 1))
                st = pr = v
        out.append((st, pr + 1))
    return out


def make_zhunt_worker(zhunt_bin, scripts_dir, win, mn, mx, threshold):
    """Замыкание-воркер для одного куска (через threads, без pickling)."""

    def zhunt_chunk(task):
        chrom, gstart, owned_end, sub = task
        acgt = sum(sub.count(b) for b in "ACGT")
        if acgt == 0:
            return chrom, []
        if acgt == len(sub):  # нет N -> прямая координата
            seq_clean = sub
            gmap = lambda i: gstart + i
        else:  # есть N -> маппинг через индексы
            keep = [j for j, b in enumerate(sub) if b in "ACGT"]
            seq_clean = "".join(sub[j] for j in keep)
            gmap = lambda i: gstart + keep[i]
        with tempfile.NamedTemporaryFile(
            "w", suffix=".seq", dir=str(scripts_dir), delete=False
        ) as tf:
            tf.write(seq_clean)
            tmp = tf.name
        try:
            subprocess.run(
                [str(zhunt_bin), str(win), str(mn), str(mx), tmp],
                check=True,
                stdin=subprocess.DEVNULL,
                stdout=subprocess.DEVNULL,
            )
            out = subprocess.run(
                [
                    "awk",
                    "-v",
                    f"t={threshold}",
                    "NR>1 && $3>t {print NR-2}",
                    tmp + ".Z-SCORE",
                ],
                check=True,
                capture_output=True,
                text=True,
            ).stdout
        finally:
            for f in (tmp, tmp + ".Z-SCORE"):
                if os.path.exists(f):
                    os.remove(f)
        flagged = []
        for x in out.split():
            p = gmap(int(x))
            if p < owned_end:  # выбросить хвост из зоны OVERLAP
                flagged.append(p)
        return chrom, flagged

    return zhunt_chunk


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--fasta", type=Path, default=BASE_DIR / "data" / "genome.fa")
    p.add_argument(
        "--g4-bed", type=Path, default=BASE_DIR / "results" / "bed" / "g4.bed"
    )
    p.add_argument(
        "--zhunt-bed", type=Path, default=BASE_DIR / "results" / "bed" / "zhunt.bed"
    )
    p.add_argument("--scripts-dir", type=Path, default=BASE_DIR / "scripts")
    p.add_argument(
        "--chroms",
        nargs="*",
        default=None,
        help="имена хромосом; не указано или 'all' = весь геном",
    )
    p.add_argument(
        "--max-bp", default="None", help="bp на хромосому; 'None'/'all' = вся длина"
    )
    p.add_argument("--only", choices=["both", "g4", "zhunt"], default="both")
    p.add_argument("--zhunt-win", type=int, default=12)
    p.add_argument("--zhunt-min", type=int, default=8)
    p.add_argument("--zhunt-max", type=int, default=12)
    p.add_argument("--z-threshold", type=float, default=400.0)
    p.add_argument("--zhunt-chunk", type=int, default=2_000_000)
    p.add_argument("--zhunt-overlap", type=int, default=1_000)
    p.add_argument("--zhunt-proc", type=int, default=96)
    return p.parse_args()


def main():
    args = parse_args()
    from Bio import SeqIO

    chroms = None if (not args.chroms or args.chroms == ["all"]) else set(args.chroms)
    max_bp = None if str(args.max_bp).lower() in ("none", "all") else int(args.max_bp)

    args.g4_bed.parent.mkdir(parents=True, exist_ok=True)
    if not args.fasta.exists():
        raise FileNotFoundError(f"нет генома {args.fasta}")

    do_g4 = args.only in ("both", "g4")
    do_zhunt = args.only in ("both", "zhunt")

    zhunt_bin = None
    worker = None
    if do_zhunt:
        zhunt_bin = build_zhunt(args.scripts_dir)
        print("zhunt собран:", zhunt_bin)
        worker = make_zhunt_worker(
            zhunt_bin,
            args.scripts_dir,
            args.zhunt_win,
            args.zhunt_min,
            args.zhunt_max,
            args.z_threshold,
        )

    print("сканирую заголовки генома для подсчёта хромосом...")
    target_ids = list_chrom_ids(args.fasta, chroms)
    total_chrom = len(target_ids)
    print(f"будет обработано хромосом: {total_chrom}")

    g4f = open(args.g4_bed, "w") if do_g4 else None
    zf = open(args.zhunt_bed, "w") if do_zhunt else None
    pool = ThreadPoolExecutor(max_workers=args.zhunt_proc) if do_zhunt else None

    n_g4 = n_zhunt = 0
    idx = 0
    try:
        for rec in SeqIO.parse(str(args.fasta), "fasta"):
            if chroms is not None and rec.id not in chroms:
                continue
            idx += 1
            s = str(rec.seq).upper()
            if max_bp is not None:
                s = s[:max_bp]
            L = len(s)
            print(f"[{idx}/{total_chrom}] хромосома {rec.id} (len={L})")

            if do_g4:
                for c, st, en, strand in find_g4(s, rec.id):
                    g4f.write(f"{c}\t{st}\t{en}\tG4\t{en - st}\t{strand}\n")
                    n_g4 += 1
                g4f.flush()

            if do_zhunt:
                tasks = []
                for st in range(0, L, args.zhunt_chunk):
                    owned_end = min(st + args.zhunt_chunk, L)
                    en = min(st + args.zhunt_chunk + args.zhunt_overlap, L)
                    tasks.append((rec.id, st, owned_end, s[st:en]))
                positions = []
                for _chrom, flagged in tqdm(pool.map(worker, tasks),
                                            total=len(tasks),
                                            desc=f"zhunt {rec.id} [{idx}/{total_chrom}]"):
                    positions.extend(flagged)
                for st, en in merge_positions(positions):
                    zf.write(f"{rec.id}\t{st}\t{en}\tZHUNT\t{en - st}\t.\n")
                    n_zhunt += 1
                zf.flush()

            print(f"  -> {rec.id}: G4={n_g4} | ZHUNT={n_zhunt}")
    finally:
        if g4f:
            g4f.close()
        if zf:
            zf.close()
        if pool:
            pool.shutdown()

    if do_g4:
        print("готово ->", args.g4_bed, "| G4-интервалов:", n_g4)
    if do_zhunt:
        print("готово ->", args.zhunt_bed, "| zhunt-интервалов:", n_zhunt)


if __name__ == "__main__":
    main()
