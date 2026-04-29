from pathlib import Path


def next_available_filename(path: str | Path) -> Path:
    path = Path(path)
    parent = path.parent
    stem = path.stem
    suffix = path.suffix

    if not (parent / f"{stem}{suffix}").exists():
        return parent / f"{stem}{suffix}"

    existing = list(parent.glob(f"{stem}*{suffix}"))
    numbers = [1]

    for file_path in existing:
        name = file_path.stem

        if name.startswith(stem + "(") and name.endswith(")"):
            try:
                number = int(name[len(stem) + 1 : -1])
                numbers.append(number)
            except ValueError:
                continue

    next_number = max(numbers) + 1
    return parent / f"{stem}({next_number}){suffix}"