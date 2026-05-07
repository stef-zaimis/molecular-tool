TXT_OUTPUT_BASENAME = "DMCs_output.txt"
XLSX_OUTPUT_BASENAME = "comparison_output.xlsx"

COLORS = {
    "A": "90EE90",
    "C": "87CEFA",
    "T": "FF7F7F",
    "G": "FFA500",

    # Common ambiguity states
    "R": "FFD27F",  # A/G
    "Y": "D8BFD8",  # C/T
    "W": "FFF2B2",  # A/T
    "M": "BFEFFF",  # A/C
    "S": "D9EAD3",  # C/G
    "K": "FCE5CD",  # G/T

    # Three-base ambiguity states
    "B": "EADCF8",  # C/G/T
    "D": "F4CCCC",  # A/G/T
    "H": "D0E0E3",  # A/C/T
    "V": "D9D2E9",  # A/C/G

    # Any base / missing / gap
    "N": "C0C0C0",
    "-": "E0E0E0",
    "?": "E0E0E0",
}

IUPAC = {
    # Strict DNA bases
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},

    # Two-base ambiguity codes
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "W": {"A", "T"},
    "M": {"A", "C"},
    "S": {"C", "G"},
    "K": {"G", "T"},

    # Three-base ambiguity codes
    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},

    # Any base
    "N": {"A", "C", "G", "T"},

    # Empty / missing states
    "-": set(),
    "?": set(),
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

BD_AMBIGUOUS_FRACTION_LIMIT = 0.20