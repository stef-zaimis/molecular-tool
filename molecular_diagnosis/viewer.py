import tkinter as tk
import tkinter.font as tkfont

from molecular_diagnosis.constants import COLORS


def open_fasta_viewer(root: tk.Tk, sequences: dict[str, str]) -> None:
    viewer = tk.Toplevel(root)
    viewer.title("FASTA Viewer")
    viewer.geometry("1200x500")

    font = tkfont.nametofont("TkDefaultFont")
    bold_font = font.copy()
    bold_font.configure(weight="bold")

    name_width = max(
        font.measure("Sequence name"),
        max(font.measure(header) for header in sequences),
    ) + 6

    base_width = font.measure("W") + 6
    row_height = 24
    seq_length = max(len(sequence) for sequence in sequences.values())

    position_digits = len(str(seq_length))
    digit_spacing = font.metrics("linespace") * 0.75
    header_height = int(position_digits * digit_spacing + 8)

    outer = tk.Frame(viewer)
    outer.pack(fill="both", expand=True)

    corner_canvas = tk.Canvas(
        outer,
        width=name_width,
        height=header_height,
        highlightthickness=0,
    )
    top_canvas = tk.Canvas(outer, height=header_height, highlightthickness=0)
    name_canvas = tk.Canvas(outer, width=name_width, highlightthickness=0)
    seq_canvas = tk.Canvas(outer, highlightthickness=0)

    y_scroll = tk.Scrollbar(outer, orient="vertical")
    x_scroll = tk.Scrollbar(outer, orient="horizontal")

    corner_canvas.grid(row=0, column=0, sticky="nsew")
    top_canvas.grid(row=0, column=1, sticky="ew")
    name_canvas.grid(row=1, column=0, sticky="ns")
    seq_canvas.grid(row=1, column=1, sticky="nsew")
    y_scroll.grid(row=1, column=2, sticky="ns")
    x_scroll.grid(row=2, column=1, sticky="ew")

    outer.grid_rowconfigure(1, weight=1)
    outer.grid_columnconfigure(1, weight=1)

    def yview(*args: object) -> None:
        name_canvas.yview(*args)
        seq_canvas.yview(*args)

    def xview(*args: object) -> None:
        top_canvas.xview(*args)
        seq_canvas.xview(*args)

    y_scroll.config(command=yview)
    x_scroll.config(command=xview)

    name_canvas.config(yscrollcommand=y_scroll.set)
    seq_canvas.config(xscrollcommand=x_scroll.set, yscrollcommand=y_scroll.set)

    corner_canvas.create_rectangle(
        0,
        0,
        name_width,
        header_height,
        fill="white",
        outline="gray",
    )
    corner_canvas.create_text(
        3,
        header_height / 2,
        text="Sequence name",
        anchor="w",
        font=bold_font,
    )

    for pos in range(seq_length):
        x1 = pos * base_width
        x2 = x1 + base_width
        label = f"{pos + 1:0{position_digits}d}"

        top_canvas.create_rectangle(
            x1,
            0,
            x2,
            header_height,
            fill="white",
            outline="gray",
        )

        center_x = (x1 + x2) / 2
        start_y = (header_height - ((position_digits - 1) * digit_spacing)) / 2
        seen_nonzero = False

        for digit_index, digit in enumerate(label):
            y = start_y + digit_index * digit_spacing

            if digit != "0":
                seen_nonzero = True

            digit_fill = "light gray" if digit == "0" and not seen_nonzero else "black"

            top_canvas.create_text(
                center_x,
                y,
                text=digit,
                anchor="center",
                font=font,
                fill=digit_fill,
            )

    row_text_items: dict[int, list[tuple[tk.Canvas, int]]] = {}
    hovered_row: int | None = None

    for row_index, (header, sequence) in enumerate(sequences.items()):
        y1 = row_index * row_height
        y2 = y1 + row_height

        name_canvas.create_rectangle(
            0,
            y1,
            name_width,
            y2,
            fill="white",
            outline="gray",
        )
        name_text = name_canvas.create_text(
            3,
            (y1 + y2) / 2,
            text=header,
            anchor="w",
            font=font,
        )

        padded_sequence = sequence.ljust(seq_length, "-")
        row_text_items[row_index] = [(name_canvas, name_text)]

        for pos, base in enumerate(padded_sequence):
            x1 = pos * base_width
            x2 = x1 + base_width
            fill = "#" + COLORS.get(base, "FFFFFF")

            seq_canvas.create_rectangle(
                x1,
                y1,
                x2,
                y2,
                fill=fill,
                outline="gray",
            )
            base_text = seq_canvas.create_text(
                (x1 + x2) / 2,
                (y1 + y2) / 2,
                text=base,
                anchor="center",
                font=font,
            )

            row_text_items[row_index].append((seq_canvas, base_text))

    total_height = len(sequences) * row_height
    total_width = seq_length * base_width

    top_canvas.config(scrollregion=(0, 0, total_width, header_height))
    name_canvas.config(scrollregion=(0, 0, name_width, total_height))
    seq_canvas.config(scrollregion=(0, 0, total_width, total_height))

    def set_row_font(row_index: int, selected: bool) -> None:
        if row_index not in row_text_items:
            return

        chosen_font = bold_font if selected else font

        for canvas, item in row_text_items[row_index]:
            canvas.itemconfigure(item, font=chosen_font)

    def clear_hover() -> None:
        nonlocal hovered_row

        if hovered_row is not None:
            set_row_font(hovered_row, False)
            hovered_row = None

    def row_from_event(canvas: tk.Canvas, event: tk.Event) -> int | None:
        y = canvas.canvasy(event.y)
        row_index = int(y // row_height)

        if row_index < 0 or row_index >= len(sequences):
            return None

        return row_index

    def on_row_motion(canvas: tk.Canvas, event: tk.Event) -> str:
        nonlocal hovered_row

        row_index = row_from_event(canvas, event)

        if row_index is None:
            clear_hover()
            return "break"

        if row_index != hovered_row:
            clear_hover()
            set_row_font(row_index, True)
            hovered_row = row_index

        return "break"

    def scroll_vertical(units: int) -> None:
        name_canvas.yview_scroll(units, "units")
        seq_canvas.yview_scroll(units, "units")
        clear_hover()

    def scroll_horizontal(units: int) -> None:
        top_canvas.xview_scroll(units, "units")
        seq_canvas.xview_scroll(units, "units")
        clear_hover()

    def on_mousewheel(event: tk.Event) -> str:
        if event.state & 0x0004:
            scroll_horizontal(int(-1 * (event.delta / 120)))
        else:
            scroll_vertical(int(-1 * (event.delta / 120)))

        return "break"

    def on_linux_scroll_up(event: tk.Event) -> str:
        if event.state & 0x0004:
            scroll_horizontal(-3)
        else:
            scroll_vertical(-3)

        return "break"

    def on_linux_scroll_down(event: tk.Event) -> str:
        if event.state & 0x0004:
            scroll_horizontal(3)
        else:
            scroll_vertical(3)

        return "break"

    name_canvas.bind("<Motion>", lambda event: on_row_motion(name_canvas, event))
    seq_canvas.bind("<Motion>", lambda event: on_row_motion(seq_canvas, event))

    for canvas in (corner_canvas, top_canvas, name_canvas, seq_canvas):
        canvas.bind("<MouseWheel>", on_mousewheel)
        canvas.bind("<Button-4>", on_linux_scroll_up)
        canvas.bind("<Button-5>", on_linux_scroll_down)
        canvas.bind("<Leave>", lambda _event: clear_hover())

    viewer.bind("<Leave>", lambda _event: clear_hover())