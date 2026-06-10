---
name: Eledone Individual Pipeline
overview: "Индивидуальная часть для Eledone cirrhosa (GCA_964016885.1) с готовой аннотацией Ensembl: скачивание данных, поиск 10 эпигенетических семейств через hmmer, аннотация генома вторичными структурами (квадруплексы, zhunt, ZDNABERT на подмножестве), разметка участков и таблицы распределения с фоном, оформление README на GitHub."
todos:
  - id: setup
    content: Создать структуру репозитория и conda-окружение (hmmer, bedtools, seqkit, samtools, agat, python, pytorch/transformers)
    status: pending
  - id: download
    content: "Скачать из Ensembl: genome/unmasked.fa.gz, geneset/2024_11/genes.gff3.gz и pep.fa.gz; сделать samtools faidx"
    status: pending
  - id: stats
    content: Посчитать статистику генома (seqkit stats, N50, GC, число генов) для слайда/README
    status: pending
  - id: hmmer
    content: "(основа: Тетрадка 2 HMMER) Скачать 10 Pfam HMM (PF00145, PF01429, PF12851, PF00856, PF02373, PF00850, PF01853, PF00439, PF00385, PF00125), прогнать hmmsearch по pep.fa, собрать таблицу семейство-ген-координаты"
    status: pending
  - id: g4
    content: "(основа: Тетрадка 1) Python-скрипт: G-квадруплексы regex на обоих стрендах (uppercase) -> results/bed/g4.bed"
    status: pending
  - id: zhunt
    content: "(основа: Тетрадка 1) Собрать и запустить zhunt по хромосомам, фильтр z-score>400 -> results/bed/zhunt.bed"
    status: pending
  - id: zdnabert
    content: Запустить ZDNABERT на CPU на 1-2 самых длинных хромосомах кусками -> results/bed/zdnabert.bed (пилот)
    status: pending
  - id: features
    content: "Из GFF построить BED участков: exons, introns, promoters (1000bp up TSS), downstream (200bp), intergenic"
    status: pending
  - id: intersect
    content: "bedtools intersect: Таблица 1 (число/доля структур) и Таблица 2 (число/доля участков со структурой) + сравнение с фоном (bedtools shuffle)"
    status: pending
  - id: report
    content: "Оформить README.md: описание организма, таблицы, код; закоммитить геном, bed-файлы и scripts на GitHub"
    status: pending
isProject: false
---

# Индивидуальная часть: Eledone cirrhosa (GCA_964016885.1)

Организм уже подтверждён, аннотация Ensembl (geneset 2024_11) есть. Геном ~3 Гб, 27 хромосом (`1`..`26`,`Z`), 12 965 белок-кодирующих генов, GC 36.5%. Железо: только CPU -> quadruplex/zhunt на всём геноме, ZDNABERT на 1-2 хромосомах (пилот).

## 0. Структура репозитория и окружение
```
data/         # геном, аннотация, протеом, Pfam HMM
results/bed/  # bed-файлы структур и участков
results/tables/
scripts/
README.md
```
Окружение (conda/mamba): `hmmer bedtools seqkit samtools agat gffread python=3.11 pandas biopython`; для ZDNABERT — `pytorch transformers` (CPU-сборка).

## 1. Скачивание данных (wget-стиль как в Тетрадке 1, но URL -> Ensembl)
Метод тот же, что в Тетрадке 1 (`wget` + `gunzip`); меняется только источник: у Eledone аннотации на NCBI НЕТ, поэтому геном+GFF+протеом берём из Ensembl (одни и те же имена хромосом `1..26,Z`, координаты согласованы).
```bash
BASE=https://ftp.ebi.ac.uk/pub/ensemblorganisms/Eledone_cirrhosa/GCA_964016885.1
wget $BASE/genome/unmasked.fa.gz              # геном (unmasked: поиск структур чувствителен к регистру)
wget $BASE/ensembl/geneset/2024_11/genes.gff3.gz   # аннотация
wget $BASE/ensembl/geneset/2024_11/pep.fa.gz       # протеом
gunzip *.gz
samtools faidx unmasked.fa                    # -> .fai = chrom sizes для bedtools
```
Нельзя "как в тетрадке" взять GFF/feature_table с NCBI — у Eledone там 0 генов. Геном с NCBI тоже не мешать с Ensembl GFF (имена хромосом не совпадут: NCBI `OZ…` vs Ensembl `1,2,…`).

## 2. Статистика генома (для слайда/README)
`seqkit stats -a genome.fa`; число генов из GFF (`awk '$3=="gene"'`). Зафиксировать: длину, N50, GC, число генов.

## 3. Гены эпигенетики: 10 семейств через hmmer по протеому
Основа: референс-тетрадка HMMER (см. раздел "Референс-тетрадки").
Выбранные 10 семейств (Pfam) из листа selected, покрывают ДНК-метил., модиф. гистонов и гистон-подобные:
- DNA C5-метилтрансфераза (DNMT1/3) — `PF00145`
- Methyl-CpG-binding (MBD/MeCP2) — `PF01429`
- TET-диоксигеназа (деметилирование ДНК) — `PF12851`
- SET-домен гистоновых метилтрансфераз (KMT/EZH) — `PF00856`
- JmjC деметилазы гистонов (KDM) — `PF02373`
- Гистондеацетилазы (HDAC) — `PF00850`
- MYST-ацетилтрансферазы (HAT) — `PF01853`
- Бромодомен (ридер ацетил-лизина) — `PF00439`
- Хромодомен (ридер метил-лизина; HP1/Polycomb) — `PF00385`
- Гистоновый фолд (гистон-подобные/варианты) — `PF00125`

Скачать профили (`https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/PFxxxxx?annotation=hmm` или из `Pfam-A.hmm`), затем:
```bash
hmmsearch --cut_ga --domtblout PF00145.tbl PF00145.hmm pep.fa > /dev/null
```
Спарсить hit-ы (e-value/score), сопоставить protein_id -> ген/координаты через GFF -> таблица `семейство | ген | координаты`.

## 4. Вторичные структуры ДНК -> 3 BED
Основа для G4 и zhunt: референс-тетрадка вторичных структур (см. раздел "Референс-тетрадки"). ZDNABERT пишется отдельно.
- G-квадруплексы (Python regex, оба стренда, верхний регистр):
  плюс: `(G{3,5}[ATGC]{1,7}){3,}G{3,5}`; минус: тот же по C (`(C{3,5}...){3,}C{3,5}`) -> `results/bed/g4.bed`.
- zhunt: собрать C-код, прогнать по хромосомам, фильтр `z-score > 400` -> `results/bed/zhunt.bed`.
- ZDNABERT (CPU): на 1-2 самых длинных хромосомах, нарезка на окна/k-меры кусками (избежать OOM) -> `results/bed/zdnabert.bed` (пилот; явно отметить в README, что подмножество).

## 5. Разметка участков генома из GFF -> BED
Через `agat`/собственный скрипт + chrom sizes:
- `exons.bed`, `introns.bed`
- `promoters.bed` = 1000 bp upstream от TSS (учесть стренд)
- `downstream.bed` = 200 bp после конца гена
- `intergenic.bed` = `bedtools complement` от объединения генов

## 6. Пересечения -> 2 таблицы + фон
- Таблица 1 (число и доля структур по участкам): `bedtools intersect -a structure -b feature` -> подсчёт; доля = структур в участке / всего структур.
- Таблица 2 (число и доля участков, содержащих структуру): `bedtools intersect -c` по участкам, доля = участков с >=1 структурой / всего участков.
- Фон: `bedtools shuffle` (рандомные регионы той же длины по геному, `-g chrom.sizes`) или сравнение с геном-wide ожидаемой долей; сравнить с Таблицей 1.

## 7. Отчёт (README.md на GitHub)
- Краткое описание Eledone cirrhosa.
- В репозитории: файл генома + bed-файлы результатов.
- Таблица распределения структур (Таблицы 1 и 2 + фон).
- Таблица генов эпигенетики (10 семейств, п.3).
- Используемый код (`scripts/`).

## Референс-тетрадки (Colab)
Готовые тетрадки, на которых строятся отдельные шаги:

- Тетрадка HMMER — https://colab.research.google.com/drive/1fijxgbCtMDlCZcpSPCFm0C4eDhLqxBFy
  Покрывает: шаг 3 (todo `hmmer`). Установка HMMER, скачивание `Pfam-A.hmm` (Pfam 35.0), `hmmpress`/`hmmfetch` профиля, `hmmsearch` по протеому.
  Адаптация под Eledone: искать по `pep.fa` Eledone (а не по протеому человека); вытащить 10 нужных PF с версиями (`hmmfetch PFxxxxx.N`); добавить `--domtblout` + `--cut_ga` для машинной таблицы.

- Тетрадка вторичных структур — https://colab.research.google.com/drive/1hHs0m3989O-PF9XFHXwu0RILa_k7TZId
  Покрывает: шаг 4 частично (todo `g4`, `zhunt`). Regex-квадруплексы на обоих стрендах (Biopython) + zhunt (`zhunt3-alan.c`, чтение `.Z-SCORE`, фильтр).
  Адаптация под Eledone: источник данных Ensembl/Eledone (в тетрадке — Octopus bimaculoides с NCBI); привести паттерн к `G{3,5}` (в тетрадке `G{3,}`); порог z-score `>400` (в тетрадке `>=500`); фикс координат минус-цепи (лишний пересчёт `seq_length - end`).
  НЕ покрывает: ZDNABERT (todo `zdnabert`), разметку участков из GFF (todo `features`), bedtools-таблицы и фон (todo `intersect`).

Покрытие по todo: `hmmer` <- Тетрадка 2; `g4`, `zhunt` <- Тетрадка 1; `download` (частично, шаблон wget) <- Тетрадка 1. Остальное (`zdnabert`, `features`, `intersect`, `report`, `setup`, `stats`) — пишется отдельно.

## Заметки по рискам
- ZDNABERT на 3 Гб без GPU нереально целиком -> только подмножество, честно указать параметры и охват.
- Промоторы/downstream обрезать по границам хромосом (`bedtools slop`/`flank` + `-g`).
- Для слайда число генов = 12 965 (Ensembl 2024_11); идентификатор = `GCA_964016885.1`.