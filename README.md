# Eledone cirrhosa — эпигенетика и вторичные структуры ДНК

Поиск генов эпигенетической машины и вторичных структур ДНК (Z-ДНК, G-квадруплексы) в геноме
курчавого осьминога *Eledone cirrhosa*.

## Организм

- **Вид:** курчавый осьминог (*Eledone cirrhosa* Lamarck, 1798; англ. curled /
  horned octopus). Семейство Octopodidae, подсемейство Eledoninae.
- **Среда обитания:** Северо-восточная Атлантика и Средиземное море; бентосный
  вид на мягких грунтах, глубины 50–700 м; ~8–18 °C (оптимум 10–14 °C).
- **Размер генома:** 2.996 Gb (2 995 825 281 bp), GC 36.5%;
  6 682 контига → 1 795 скаффолдов. PacBio HiFi 32x + Arima2 Hi-C
  (Hifiasm + purge_dups + YaHS + ручная курация TreeVal). Митогеном 16 147 bp.
- **N50:** Scaffold N50 = 378.8 Mb (топ-1 в отряде Octopoda; L50 = 4);
  Contig N50 = 0.82 Mb.
- **Гены:** 12 965 белок-кодирующих генов (21 014 транскриптов;
  Ensembl geneset 2024_11).
- **Публикаций в PubMed по виду:** 44.

### Эпигенетический контекст

В созревающих сперматозоидах гистоны полностью заменяются на маленький
цистеин-богатый протамин, сшивающий ДНК дисульфидными связями
(Gimenez-Bonafé 2002) — классическое для головоногих переписывание хроматина.
У coleoid-цефалопод общее CpG-метилирование низкое (<10% CpG), 5mC
сосредоточено в телах активно работающих генов; метки
H3K4me3 / H3K9me3 / H3K27me3 консервативные. Сборка получена с Arima2
Hi-C-данными, что позволяет анализировать 3D-структуру генома.

## Данные

**Геном** — `https://ftp.ebi.ac.uk/pub/ensemblorganisms/Eledone_cirrhosa/GCA_964016885.1/genome/`
Геном (~3 ГБ) и веса моделей в репозиторий не коммитятся — качаются отдельно:

```bash
mkdir -p data && cd data
wget https://ftp.ebi.ac.uk/pub/ensemblorganisms/Eledone_cirrhosa/GCA_964016885.1/genome/unmasked.fa.gz
wget https://ftp.ebi.ac.uk/pub/ensemblorganisms/Eledone_cirrhosa/GCA_964016885.1/ensembl/geneset/2024_11/genes.gff3.gz
wget https://ftp.ebi.ac.uk/pub/ensemblorganisms/Eledone_cirrhosa/GCA_964016885.1/ensembl/geneset/2024_11/pep.fa.gz
gunzip *.gz
samtools faidx unmasked.fa
```