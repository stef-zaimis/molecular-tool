TXT_OUTPUT_BASENAME = "DMCs_output.txt"
XLSX_OUTPUT_BASENAME = "comparison_output.xlsx"

COLORS = {
    "A": "90EE90",
    "C": "87CEFA",
    "T": "FF7F7F",
    "G": "FFA500",
    "N": "C0C0C0",
    "R": "FFD27F",
    "Y": "D8BFD8",
    "W": "FFF2B2",
    "M": "BFEFFF",
    "-": "E0E0E0",
}

IUPAC = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "W": {"A", "T"},
    "M": {"A", "C"},
    "N": {"A", "C", "G", "T"},
    "-": set(),
}

STRICT_BASES = {"A", "C", "G", "T"}

# Focal punishment / anomaly scoring thresholds.
# Practical boundary convention:
#   empty_fraction < 1/3      -> polymorphism
#   1/3 <= empty_fraction < 2/3 -> balancing
#   empty_fraction >= 2/3     -> prolongation/insertion logic
POLYMORPHISM_EMPTY_LIMIT = 1 / 3
BALANCING_EMPTY_LIMIT = 2 / 3

POLYMORPHISM_EMPTY_WEIGHT = 0.20
BALANCING_EMPTY_WEIGHT = 0.10
PROLONGATION_WEIGHT = 0.10