"""Воркер для параллельного инференса ZDNABERT.

Используется тетрадкой notebooks/zdnabert_eledone.ipynb через multiprocessing (spawn):
функции должны жить в импортируемом модуле, иначе spawn-процессы их не подхватят.

Модель грузится один раз на процесс в init_worker() и кладётся в глобал,
чтобы не перезагружать её на каждый чанк.
"""
import os
import numpy as np

_MODEL = None
_TOK = None


def init_worker(model_dir, nthreads):
    """Запускается один раз в каждом процессе пула (initializer)."""
    import torch
    torch.set_num_threads(nthreads)
    from transformers import BertTokenizer, BertForTokenClassification
    global _MODEL, _TOK
    _TOK = BertTokenizer.from_pretrained(model_dir)
    _MODEL = BertForTokenClassification.from_pretrained(model_dir).eval()


def seq2kmer(seq, k=6):
    return [seq[x:x + k] for x in range(len(seq) + 1 - k)]


def split_seq(seq, length=512, pad=16):
    res = []
    for st in range(0, len(seq), length - pad):
        res.append(seq[st:min(st + length, len(seq))])
    return res


def stitch_np_seq(np_seqs, pad=16):
    res = np.array([])
    for s in np_seqs:
        res = res[:-pad] if len(res) else res
        res = np.concatenate([res, s])
    return res


def predict_track(sub_seq, batch_size=64):
    """P(Z-DNA) на каждый нуклеотид куска. Батч с attention_mask -> идентично поштучному."""
    import torch
    pieces = split_seq(seq2kmer(sub_seq.upper(), 6))
    encoded = [_TOK.encode(" ".join(p), add_special_tokens=False) for p in pieces]
    preds = []
    with torch.no_grad():
        for i in range(0, len(encoded), batch_size):
            batch = encoded[i:i + batch_size]
            maxlen = max(len(x) for x in batch)
            input_ids = torch.zeros(len(batch), maxlen, dtype=torch.long)
            attn = torch.zeros(len(batch), maxlen, dtype=torch.long)
            for j, x in enumerate(batch):
                input_ids[j, :len(x)] = torch.tensor(x, dtype=torch.long)
                attn[j, :len(x)] = 1
            probs = torch.softmax(_MODEL(input_ids, attention_mask=attn)[-1], dim=-1)[:, :, 1]
            for j, x in enumerate(batch):
                preds.append(probs[j, :len(x)].numpy())
    return stitch_np_seq(preds)


def process_chunk(task):
    """task = (chrom, offset, seq_chunk, threshold, min_len, batch_size, cache_path).

    Возвращает список (chrom, start, end) в глобальных координатах хромосомы.
    Если cache_path существует — читает результат оттуда (модель не трогает).
    Иначе считает, атомарно пишет файл-кэш (даже пустой) и возвращает результат.
    """
    chrom, offset, seq_chunk, threshold, min_len, batch_size, cache_path = task

    # кэш-хит: читаем готовый результат
    if cache_path is not None and os.path.exists(cache_path):
        out = []
        with open(cache_path) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    continue
                c, s, e = line.split("\t")
                out.append((c, int(s), int(e)))
        return out

    # счёт
    import scipy.ndimage as ndi
    if len(seq_chunk) < 6:
        out = []
    else:
        track = predict_track(seq_chunk, batch_size)
        labeled, n = ndi.label(track > threshold)
        out = []
        for lab in range(1, n + 1):
            idx = np.where(labeled == lab)[0]
            if idx.shape[0] > min_len:
                out.append((chrom, offset + int(idx[0]), offset + int(idx[-1]) + 1))

    # атомарная запись (пустой файл = валидный маркер "посчитано, 0 интервалов")
    if cache_path is not None:
        tmp = cache_path + ".tmp"
        with open(tmp, "w") as fh:
            for c, s, e in out:
                fh.write(f"{c}\t{s}\t{e}\n")
        os.replace(tmp, cache_path)
    return out
