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