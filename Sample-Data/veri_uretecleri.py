"""Biyoenformatiğe Giriş dersi - sentetik veri seti üreteçleri.

Her fonksiyon kendi içinde bağımsızdır (importlar dahil) ve notebook'lara
aynen gömülür. Tohum (seed) sabit olduğundan her çalıştırmada aynı veri üretilir.
"""


def veri_uret_plazmit(klasor="veri"):
    """ornek_plazmit.gb (GenBank) ve ornek_genler.fasta dosyalarını üretir."""
    import os, random
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from Bio import SeqIO
    os.makedirs(klasor, exist_ok=True)
    rng = random.Random(42)
    stop = ["TAA", "TAG", "TGA"]
    anlamli = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"
               if a + b + c not in stop]
    engel = ["GAATTC", "GGATCC", "AAGCTT"]          # EcoRI, BamHI, HindIII

    def rastgele(n, gc=0.5):
        w = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
        while True:
            s = "".join(rng.choices("ACGT", weights=w, k=n))
            if not any(e in s for e in engel):
                return s

    def orf(n_kodon, gc=0.5):
        while True:
            govde = "".join(rng.choice(anlamli) for _ in range(n_kodon - 2))
            s = "ATG" + govde + rng.choice(stop)
            if not any(e in s for e in engel):
                return s

    parcalar, ozellikler, pos = [], [], 0

    def ekle(s, tip=None, **q):
        nonlocal pos
        bas = pos
        parcalar.append(s)
        pos += len(s)
        if tip:
            ozellikler.append((tip, bas, pos, q))
        return bas, pos

    ekle(rastgele(180))
    ekle("TTGACA" + rastgele(17) + "TATAAT", "promoter", note="sigma70 benzeri promotor (-35/-10)")
    ekle(rastgele(12))
    ekle("AGGAGG", "RBS", note="Shine-Dalgarno dizisi")
    ekle(rastgele(7))
    cds1 = orf(300)
    ekle(cds1, "CDS+", gene="ornA", product="ornek protein A")
    ekle(rastgele(60))
    ekle("GAATTC", "misc_feature", note="EcoRI kesim yeri")
    ekle(rastgele(120))
    cds2 = orf(200, gc=0.45)
    ekle(str(Seq(cds2).reverse_complement()), "CDS-", gene="ornB", product="ornek protein B")
    ekle(rastgele(80))
    ekle("GGATCC", "misc_feature", note="BamHI kesim yeri")
    ekle(rastgele(150))
    ekle("AAGCTT", "misc_feature", note="HindIII kesim yeri")
    ekle(rastgele(3000 - pos - 20))
    ekle("GAATTC", "misc_feature", note="EcoRI kesim yeri (2)")
    ekle(rastgele(3000 - pos))
    dizi = "".join(parcalar)

    kayit = SeqRecord(Seq(dizi), id="pORNEK1", name="pORNEK1",
                      description="Sentetik ornek plazmit, Biyoenformatige Giris dersi")
    kayit.annotations["molecule_type"] = "DNA"
    kayit.annotations["topology"] = "circular"
    kayit.annotations["organism"] = "sentetik yapi"
    kayit.annotations["data_file_division"] = "SYN"
    kayit.annotations["date"] = "01-SEP-2026"
    kayit.features.append(SeqFeature(FeatureLocation(0, len(dizi), strand=1), type="source",
                                     qualifiers={"organism": ["sentetik yapi"], "mol_type": ["genomic DNA"]}))
    for tip, b, s, q in ozellikler:
        if tip.startswith("CDS"):
            strand = 1 if tip == "CDS+" else -1
            f = SeqFeature(FeatureLocation(b, s, strand=strand), type="CDS",
                           qualifiers={k: [v] for k, v in q.items()})
            f.qualifiers["translation"] = [str(f.extract(kayit.seq).translate(to_stop=True))]
            f.qualifiers["codon_start"] = ["1"]
            f.qualifiers["transl_table"] = ["11"]
        else:
            f = SeqFeature(FeatureLocation(b, s, strand=1), type=tip,
                           qualifiers={k: [v] for k, v in q.items()})
        kayit.features.append(f)
    SeqIO.write(kayit, os.path.join(klasor, "ornek_plazmit.gb"), "genbank")

    genler = [
        SeqRecord(Seq(cds1), id="ornA", description="plazmitteki ileri yonlu gen"),
        SeqRecord(Seq(cds2), id="ornB", description="plazmitteki ters yonlu gen (5'->3')"),
        SeqRecord(Seq(orf(150, gc=0.68)), id="yuksekGC", description="GC orani yuksek sentetik gen"),
        SeqRecord(Seq(orf(150, gc=0.32)), id="dusukGC", description="GC orani dusuk sentetik gen"),
        SeqRecord(Seq(orf(80)), id="kisaGen", description="kisa sentetik gen"),
    ]
    # GC eğilimli genleri kodon seçimiyle üret
    def gc_orf(n, gc):
        agirlik = []
        for k in anlamli:
            g = sum(c in "GC" for c in k)
            agirlik.append((gc / 0.5) ** g * ((1 - gc) / 0.5) ** (3 - g))
        while True:
            s = "ATG" + "".join(rng.choices(anlamli, weights=agirlik, k=n - 2)) + "TAA"
            if not any(e in s for e in engel):
                return s
    genler[2].seq = Seq(gc_orf(150, 0.68))
    genler[3].seq = Seq(gc_orf(150, 0.32))
    SeqIO.write(genler, os.path.join(klasor, "ornek_genler.fasta"), "fasta")


def veri_uret_fastq(klasor="veri"):
    """okumalar.fastq: 2000 adet 100 bazlık sentetik Illumina benzeri okuma."""
    import os, random, math
    os.makedirs(klasor, exist_ok=True)
    rng = random.Random(7)
    genom = "".join(rng.choice("ACGT") for _ in range(5000))
    adaptor = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
    L = 100
    satirlar = []
    for i in range(2000):
        kotu = rng.random() < 0.10
        adaptorlu = rng.random() < 0.08
        if adaptorlu:
            ins = rng.randint(55, 90)
            bas = rng.randint(0, len(genom) - ins)
            gercek = genom[bas:bas + ins] + adaptor
            gercek = (gercek + "".join(rng.choice("ACGT") for _ in range(L)))[:L]
        else:
            bas = rng.randint(0, len(genom) - L)
            gercek = genom[bas:bas + L]
        dizi, kal = [], []
        for j, b in enumerate(gercek):
            ort = (38 - 13 * (j / L) ** 2) if not kotu else (22 - 12 * j / L)
            q = int(round(rng.gauss(ort, 3)))
            q = max(2, min(41, q))
            if rng.random() < 0.002 or (kotu and rng.random() < 0.01):
                b, q = "N", 2
            elif rng.random() < 10 ** (-q / 10):
                b = rng.choice([x for x in "ACGT" if x != b])
            dizi.append(b)
            kal.append(chr(q + 33))
        satirlar += [f"@okuma_{i+1:04d} ornek=S1", "".join(dizi), "+", "".join(kal)]
    with open(os.path.join(klasor, "okumalar.fastq"), "w") as f:
        f.write("\n".join(satirlar) + "\n")


def veri_uret_protein_ailesi(klasor="veri"):
    """protein_ailesi.fasta (hizasız), protein_ailesi_hizali.aln (gerçek MSA, Clustal)."""
    import os, random
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Align import MultipleSeqAlignment
    from Bio import SeqIO, AlignIO
    os.makedirs(klasor, exist_ok=True)
    rng = random.Random(2026)
    aa = "ACDEFGHIKLMNPQRSTVWY"
    frek = [8.3, 1.4, 5.5, 6.8, 3.9, 7.1, 2.3, 5.9, 5.8, 9.7, 2.4, 4.1, 4.7, 3.9, 5.4, 6.6, 5.3, 6.9, 1.1, 2.9]
    # Yer değiştirme olasılıkları BLOSUM62'den türetilir: P(b|a) ∝ p_b * 2^(s(a,b)/2)
    from Bio.Align import substitution_matrices
    B62 = substitution_matrices.load("BLOSUM62")
    gecis = {a: [0 if b == a else frek[j] * 2 ** (B62[a][b] / 2) for j, b in enumerate(aa)] for a in aa}
    L = 160
    ata = ["M"] + rng.choices(aa, weights=frek, k=L - 1)
    korunmus = set(rng.sample(range(1, L), 35)) | {0}

    def evrim(dizi, t):
        yeni = list(dizi)
        for i in range(len(yeni)):
            if yeni[i] == "-" or i in korunmus:
                continue
            if rng.random() < t:
                yeni[i] = rng.choices(aa, weights=gecis[yeni[i]])[0]
        if rng.random() < t * 2.5:                      # delesyon olayı
            u = rng.randint(1, 5)
            b = rng.randint(5, L - u - 5)
            if not any(k in korunmus for k in range(b, b + u)):
                for k in range(b, b + u):
                    yeni[k] = "-"
        return yeni

    # ((A1,A2),(A3,(A4,A5))) , ((B1,B2),B3)
    kA = evrim(ata, 0.15); kB = evrim(ata, 0.25)
    kA1 = evrim(kA, 0.08); kA2 = evrim(kA, 0.12)
    diziler = {
        "protA1": evrim(kA1, 0.05), "protA2": evrim(kA1, 0.06),
        "protA3": evrim(kA2, 0.10), "protA4": evrim(evrim(kA2, 0.05), 0.05),
        "protA5": evrim(evrim(kA2, 0.05), 0.07),
        "protB1": evrim(kB, 0.08), "protB2": evrim(kB, 0.10), "protB3": evrim(kB, 0.20),
    }
    hizali = MultipleSeqAlignment([SeqRecord(Seq("".join(v)), id=k) for k, v in diziler.items()])
    AlignIO.write(hizali, os.path.join(klasor, "protein_ailesi_hizali.aln"), "clustal")
    kayitlar = [SeqRecord(Seq("".join(v).replace("-", "")), id=k, description="sentetik protein ailesi uyesi")
                for k, v in diziler.items()]
    ilgisiz = "M" + "".join(rng.choices(aa, weights=frek, k=149))
    kayitlar.append(SeqRecord(Seq(ilgisiz), id="ilgisiz1", description="aileyle akraba olmayan rastgele protein"))
    SeqIO.write(kayitlar, os.path.join(klasor, "protein_ailesi.fasta"), "fasta")


def veri_uret_blast(klasor="veri"):
    """blast_sorgu.fasta, blast_veritabani.fasta, blast_gercek_homologlar.csv"""
    import os, random, csv
    os.makedirs(klasor, exist_ok=True)
    rng = random.Random(11)
    rs = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    sorgu = rs(300)

    def mutasyon(s, oran, indel=0.0):
        out = []
        for b in s:
            r = rng.random()
            if r < indel / 2:
                continue                                  # delesyon
            if r < indel:
                out.append(rng.choice("ACGT"))            # insersiyon
            out.append(rng.choice([x for x in "ACGT" if x != b]) if rng.random() < oran else b)
        return "".join(out)

    homolog = {5: (0.05, 0.0), 23: (0.10, 0.0), 47: (0.20, 0.0), 81: (0.30, 0.0),
               102: (0.12, 0.03), 133: (0.45, 0.0)}
    kayit, gercek = [], []
    for i in range(150):
        uzun = rng.randint(400, 800)
        s = rs(uzun)
        if i in homolog:
            oran, indel = homolog[i]
            bas_q, bit_q = (0, 300) if i != 47 else (50, 250)
            parca = mutasyon(sorgu[bas_q:bit_q], oran, indel)
            yer = rng.randint(20, uzun - len(parca) - 20)
            s = s[:yer] + parca + s[yer + len(parca):]
            gercek.append((f"db_{i:03d}", oran, indel, bas_q, bit_q, yer))
        kayit.append((f"db_{i:03d}", s))
    with open(os.path.join(klasor, "blast_sorgu.fasta"), "w") as f:
        f.write(">sorgu1 sentetik sorgu dizisi\n" + "\n".join(sorgu[i:i+60] for i in range(0, 300, 60)) + "\n")
    with open(os.path.join(klasor, "blast_veritabani.fasta"), "w") as f:
        for k, s in kayit:
            f.write(f">{k}\n" + "\n".join(s[i:i+60] for i in range(0, len(s), 60)) + "\n")
    with open(os.path.join(klasor, "blast_gercek_homologlar.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["id", "mutasyon_orani", "indel_orani", "sorgu_bas", "sorgu_bit", "db_konum"])
        w.writerows(gercek)


def veri_uret_filogeni(klasor="veri"):
    """filogeni_hizali.fasta (JC69 ile ağaç boyunca simüle), gercek_agac.nwk"""
    import os, random, math
    os.makedirs(klasor, exist_ok=True)
    rng = random.Random(5)
    newick = ("((((Tur_1:0.04,Tur_2:0.04):0.08,Tur_3:0.12):0.10,(Tur_4:0.09,Tur_5:0.09):0.13):0.08,"
              "((Tur_6:0.07,Tur_7:0.07):0.11,Tur_8:0.18):0.12);")
    # küçük bir newick ayrıştırıcı
    def ayristir(s):
        s = s.strip().rstrip(";")
        def dugum(i):
            if s[i] == "(":
                cocuklar = []
                i += 1
                while True:
                    c, i = dugum(i)
                    cocuklar.append(c)
                    if s[i] == ",":
                        i += 1
                        continue
                    i += 1  # ')'
                    break
                ad = ""
            else:
                cocuklar, ad = [], ""
                while s[i] not in ":,)":
                    ad += s[i]; i += 1
            uz = 0.0
            if i < len(s) and s[i] == ":":
                i += 1; j = i
                while i < len(s) and s[i] not in ",)":
                    i += 1
                uz = float(s[j:i])
            return {"ad": ad, "uz": uz, "c": cocuklar}, i
        return dugum(0)[0]
    kok = ayristir(newick)
    L = 800
    kok_dizi = [rng.choice("ACGT") for _ in range(L)]
    sonuc = {}
    def in_(d, dizi):
        p = 0.75 * (1 - math.exp(-4 * d["uz"] / 3))
        yeni = [rng.choice([x for x in "ACGT" if x != b]) if rng.random() < p else b for b in dizi]
        if not d["c"]:
            sonuc[d["ad"]] = "".join(yeni)
        for c in d["c"]:
            in_(c, yeni)
    in_(kok, kok_dizi)
    with open(os.path.join(klasor, "filogeni_hizali.fasta"), "w") as f:
        for k in sorted(sonuc):
            f.write(f">{k}\n" + "\n".join(sonuc[k][i:i+80] for i in range(0, L, 80)) + "\n")
    with open(os.path.join(klasor, "gercek_agac.nwk"), "w") as f:
        f.write(newick + "\n")


def veri_uret_ifade(klasor="veri"):
    """ifade_sayimlar.csv (RNA-seq benzeri ham sayımlar), ornek_bilgisi.csv, gercek_de_genler.csv"""
    import os
    import numpy as np
    import pandas as pd
    os.makedirs(klasor, exist_ok=True)
    rng = np.random.default_rng(3)
    G, n = 2000, 12
    genler = [f"GEN{i:04d}" for i in range(1, G + 1)]
    ornekler = [f"K{i}" for i in range(1, 7)] + [f"T{i}" for i in range(1, 7)]
    grup = ["kontrol"] * 6 + ["tedavi"] * 6
    temel = np.exp(rng.normal(4.5, 1.6, G))              # gen başına ortalama ifade
    fc = np.ones(G)
    de = rng.choice(G, 150, replace=False)
    yukari, asagi = de[:75], de[75:]
    fc[yukari] = 2 ** rng.uniform(1, 3, 75)
    fc[asagi] = 2 ** -rng.uniform(1, 3, 75)
    kutuphane = rng.uniform(0.6, 1.5, n)
    disp = 0.08
    sayim = np.zeros((G, n), dtype=int)
    for j in range(n):
        mu = temel * kutuphane[j] * (fc if grup[j] == "tedavi" else 1)
        lam = rng.gamma(1 / disp, mu * disp)
        sayim[:, j] = rng.poisson(lam)
    pd.DataFrame(sayim, index=pd.Index(genler, name="gen"), columns=ornekler).to_csv(
        os.path.join(klasor, "ifade_sayimlar.csv"))
    pd.DataFrame({"ornek": ornekler, "grup": grup,
                  "kutuphane_carpani": np.round(kutuphane, 3)}).to_csv(
        os.path.join(klasor, "ornek_bilgisi.csv"), index=False)
    pd.DataFrame({"gen": [genler[i] for i in de],
                  "yon": ["yukari"] * 75 + ["asagi"] * 75,
                  "gercek_log2FC": np.round(np.log2(fc[de]), 3)}).to_csv(
        os.path.join(klasor, "gercek_de_genler.csv"), index=False)


def veri_uret_tumor(klasor="veri"):
    """tumor_ifade.csv (200 örnek x 1000 gen, log2 ölçek), tumor_etiketler.csv, tumor_marker_genler.csv"""
    import os
    import numpy as np
    import pandas as pd
    os.makedirs(klasor, exist_ok=True)
    rng = np.random.default_rng(12)
    G = 1000
    tipler = {"A": 70, "B": 60, "C": 45, "D": 25}
    genler = [f"g{i:04d}" for i in range(1, G + 1)]
    mu = rng.normal(7, 1.5, G)
    sd = rng.uniform(0.5, 1.1, G)
    # gen modülleri: ortak gizli faktörlerle korelasyon
    faktor_yuk = np.zeros((G, 5))
    for k in range(5):
        idx = rng.choice(G, 80, replace=False)
        faktor_yuk[idx, k] = rng.uniform(0.2, 0.5, 80)
    marker = {}
    kalan = list(rng.permutation(G))
    for t in tipler:
        marker[t] = [kalan.pop() for _ in range(30)]
    ortak_CD = [kalan.pop() for _ in range(20)]            # C ve D alt tiplerinin ortak belirteçleri
    satirlar, etiket, merkez = [], [], []
    for t, n in tipler.items():
        for _ in range(n):
            x = mu + sd * rng.normal(size=G) + faktor_yuk @ rng.normal(size=5)
            idx = marker[t]
            if t in ("A", "B"):                             # belirgin alt tipler
                x[idx[:20]] += rng.uniform(0.7, 1.4, 20)
                x[idx[20:]] -= rng.uniform(0.7, 1.4, 10)
            else:                                           # C ve D birbirine benzer, farkları zayıf
                x[ortak_CD] += rng.uniform(0.7, 1.4, 20)
                x[idx[:20]] += rng.uniform(0.25, 0.6, 20)
                x[idx[20:]] -= rng.uniform(0.25, 0.6, 10)
            m = "Merkez1" if rng.random() < 0.55 else "Merkez2"
            if m == "Merkez2":
                x[:120] += 0.6                             # yığın (batch) etkisi
            satirlar.append(x); etiket.append(t); merkez.append(m)
    X = np.round(np.array(satirlar), 3)
    sira = rng.permutation(len(X))
    ornek = [f"TMR{i:03d}" for i in range(1, len(X) + 1)]
    pd.DataFrame(X[sira], index=pd.Index(ornek, name="ornek"), columns=genler).to_csv(
        os.path.join(klasor, "tumor_ifade.csv"))
    pd.DataFrame({"ornek": ornek, "alt_tip": np.array(etiket)[sira],
                  "merkez": np.array(merkez)[sira]}).to_csv(
        os.path.join(klasor, "tumor_etiketler.csv"), index=False)
    pd.DataFrame([(genler[g], t, "yukari" if i < 20 else "asagi")
                  for t, lst in marker.items() for i, g in enumerate(lst)]
                 + [(genler[g], "C+D", "yukari") for g in ortak_CD],
                 columns=["gen", "alt_tip", "yon"]).to_csv(
        os.path.join(klasor, "tumor_marker_genler.csv"), index=False)


def veri_uret_promotor(klasor="veri"):
    """promotor_veriseti.csv: 2000 adet 80 bp dizi, etiket 1=promotor, 0=promotor değil"""
    import os, random, csv
    os.makedirs(klasor, exist_ok=True)
    rng = random.Random(21)
    tata_pwm = [  # TATAAA benzeri motif için pozisyon ağırlık matrisi (A,C,G,T)
        [0.05, 0.05, 0.05, 0.85], [0.85, 0.05, 0.05, 0.05], [0.10, 0.05, 0.05, 0.80],
        [0.80, 0.05, 0.05, 0.10], [0.70, 0.05, 0.05, 0.20], [0.75, 0.05, 0.10, 0.10]]
    def bg(n, gc):
        return "".join(rng.choices("ACGT", weights=[(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2], k=n))
    def motif():
        return "".join(rng.choices("ACGT", weights=w)[0] for w in tata_pwm)
    satirlar = []
    for i in range(2000):
        pozitif = i % 2 == 0
        if pozitif:
            s = list(bg(80, rng.uniform(0.48, 0.62)))
            if rng.random() < 0.80:
                p = rng.randint(22, 34)
                s[p:p + 6] = motif()
            if rng.random() < 0.50:
                p = rng.randint(2, 15)
                s[p:p + 6] = "GGGCGG"
        else:
            s = list(bg(80, rng.uniform(0.38, 0.55)))
            if rng.random() < 0.20:
                p = rng.randint(0, 74)
                s[p:p + 6] = motif()
        satirlar.append((f"dizi_{i+1:04d}", "".join(s), int(pozitif)))
    rng.shuffle(satirlar)
    with open(os.path.join(klasor, "promotor_veriseti.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["id", "dizi", "etiket"])
        w.writerows(satirlar)


def veri_uret_tf(klasor="veri"):
    """tf_baglanma.csv: 4000 adet 101 bp dizi; etiket 1 = TF bağlanma bölgesi içeriyor"""
    import os, random, csv
    os.makedirs(klasor, exist_ok=True)
    rng = random.Random(99)
    pwm = [  # E-box çekirdekli 10 bp motif (A,C,G,T)
        [0.2, 0.3, 0.3, 0.2], [0.1, 0.6, 0.2, 0.1], [0.05, 0.9, 0.03, 0.02], [0.9, 0.03, 0.05, 0.02],
        [0.03, 0.9, 0.05, 0.02], [0.02, 0.05, 0.9, 0.03], [0.02, 0.03, 0.05, 0.9], [0.03, 0.05, 0.9, 0.02],
        [0.1, 0.2, 0.6, 0.1], [0.2, 0.3, 0.3, 0.2]]
    def motif():
        return "".join(rng.choices("ACGT", weights=w)[0] for w in pwm)
    def rs(n):
        return "".join(rng.choice("ACGT") for _ in range(n))
    satirlar = []
    for i in range(4000):
        s = list(rs(101))
        pozitif = i % 2 == 0
        if pozitif:
            p = rng.randint(0, 91)
            s[p:p + 10] = motif()
            if rng.random() < 0.3:                          # ikinci kopya
                p = rng.randint(0, 91)
                s[p:p + 10] = motif()
        else:
            if rng.random() < 0.35:                         # yanıltıcı: motifin karıştırılmış hali
                m = list(motif()); rng.shuffle(m)
                p = rng.randint(0, 91)
                s[p:p + 10] = m
        satirlar.append((f"tf_{i+1:04d}", "".join(s), int(pozitif)))
    rng.shuffle(satirlar)
    with open(os.path.join(klasor, "tf_baglanma.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["id", "dizi", "etiket"])
        w.writerows(satirlar)


def veri_uret_peptit(klasor="veri"):
    """peptitler.csv: 600 peptit; amp=1 antimikrobiyal peptit (dengesiz: 150 pozitif)"""
    import os, random, csv
    os.makedirs(klasor, exist_ok=True)
    rng = random.Random(8)
    aa = "ACDEFGHIKLMNPQRSTVWY"
    arka = [8.3, 1.4, 5.5, 6.8, 3.9, 7.1, 2.3, 5.9, 5.8, 9.7, 2.4, 4.1, 4.7, 3.9, 5.4, 6.6, 5.3, 6.9, 1.1, 2.9]
    amp = list(arka)
    for harf, kat in {"K": 2.6, "R": 2.4, "L": 1.5, "W": 2.5, "I": 1.6, "F": 1.4, "C": 1.8,
                      "D": 0.35, "E": 0.35, "Q": 0.6, "N": 0.6}.items():
        amp[aa.index(harf)] *= kat
    satirlar = []
    for i in range(600):
        pozitif = i < 150
        n = rng.randint(12, 40)
        w = amp if pozitif else arka
        if not pozitif and rng.random() < 0.15:              # zor negatifler: katyonik ama AMP değil
            w = [x * (1.8 if a in "KR" else 1) for a, x in zip(aa, arka)]
        s = "".join(rng.choices(aa, weights=w, k=n))
        satirlar.append((f"pep_{i+1:03d}", s, int(pozitif)))
    rng.shuffle(satirlar)
    with open(os.path.join(klasor, "peptitler.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["id", "dizi", "amp"])
        w.writerows(satirlar)


TUMU = [veri_uret_plazmit, veri_uret_fastq, veri_uret_protein_ailesi, veri_uret_blast,
        veri_uret_filogeni, veri_uret_ifade, veri_uret_tumor, veri_uret_promotor,
        veri_uret_tf, veri_uret_peptit]

if __name__ == "__main__":
    import sys
    hedef = sys.argv[1] if len(sys.argv) > 1 else "veri"
    for f in TUMU:
        f(hedef)
        print("tamam:", f.__name__)
