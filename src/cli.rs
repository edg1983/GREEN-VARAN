pub const VERSION: &str = "1.4";

pub const LOGO: &str = r#"
 / __)(  _ \(  __)(  __)(  ( \ ___ / )( \ / _\ (  _ \ / _\ (  ( \
( (_ \ )   / ) _)  ) _) /    /(___)\ \/ //    \ )   //    \/    /
 \___/(__\_)(____)(____)\_)__)      \__/ \_/\_/(__\_)\_/\_/\_)__)

                                              _.-~`  `~-.
                  _.--~~~---,.__          _.,;; .   -=(@'`\
               .-`              ``~~~~--~~` ';;;       ____)
            _.'            '.              ';;;;;    '`_.'
         .-~;`               `\           ' ';;;;;__.~`
       .' .'          `'.     |           /  /;''
        \/      .---'' ``)   /'-._____.--'\  \\
       _/|    (`        /  /`              `\ \__
',    `/- \   \      __/  (_                /-\-\-`
  `;'-..___)   |     `/-\-\-`
    `-.       .'
 jgs   `~~~~``
"#;

pub const STD_CHROMS: &[&str] = &[
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11",
    "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21",
    "chr22", "chrX", "chrY", "chrM",
];

pub const DESC_ANN: &str = "Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO' ";
pub const DESC_BCSQ: &str = "Local consequence annotation from BCFtools/csq, Format: Consequence|gene|transcript|biotype|strand|amino_acid_change|dna_change";
pub const DESC_CSQ: &str = "Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position";

pub const DESC_ID: &str = "GREENDB region IDs";
pub const DESC_STDTYPE: &str = "Regulatory region type as reported in GREENDB";
pub const DESC_DBSOURCE: &str = "Sources of this GREENDB annotation";
pub const DESC_GENES: &str = "Genes potentially controlled by regulatory regions according to GREENDB";
pub const DESC_CONSTRAINT: &str = "Maximum constraint value for GREENDB regions overlapping this variant";
pub const DESC_LEVEL: &str = "GREENDB prioritization level. 0:overlap, 1:0+rare, 2:1+TFBS-DNase-UCNE, 3:2+high score,4:3+high constraint";
pub const DESC_MORE_SUPPORT: &str = "GREENDB annotation";
pub const DESC_VOI: &str = "Variant affecting a region connected to a gene of interest in GREENDB";

#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub enum Impact {
    #[value(name = "HIGH")]
    High,
    #[value(name = "MODERATE")]
    Moderate,
    #[value(name = "LOW")]
    Low,
    #[value(name = "MODIFIER")]
    Modifier,
}

impl Impact {
    pub fn as_str(self) -> &'static str {
        match self {
            Impact::High => "HIGH",
            Impact::Moderate => "MODERATE",
            Impact::Low => "LOW",
            Impact::Modifier => "MODIFIER",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub enum PrioritizationStrategy {
    #[value(name = "pileup")]
    Pileup,
    #[value(name = "levels")]
    Levels,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub enum ConnectionMode {
    #[value(name = "all")]
    All,
    #[value(name = "closest")]
    Closest,
    #[value(name = "annotated")]
    Annotated,
}

impl ConnectionMode {
    pub fn as_str(self) -> &'static str {
        match self {
            ConnectionMode::All => "all",
            ConnectionMode::Closest => "closest",
            ConnectionMode::Annotated => "annotated",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub enum GenomeBuild {
    #[value(name = "GRCh37")]
    GRCh37,
    #[value(name = "GRCh38")]
    GRCh38,
}

impl GenomeBuild {
    pub fn as_str(self) -> &'static str {
        match self {
            GenomeBuild::GRCh37 => "GRCh37",
            GenomeBuild::GRCh38 => "GRCh38",
        }
    }
}

/// Resolve thread count: 0 means auto-detect available parallelism (cgroup-aware).
pub fn resolve_threads(requested: usize) -> usize {
    if requested == 0 {
        std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1)
    } else {
        requested
    }
}
