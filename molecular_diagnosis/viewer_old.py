import tkinter as tk
import tkinter.font as tkfont
from collections.abc import Mapping

from molecular_diagnosis.constants import COLORS


def open_fasta_viewer(
    root: tk.Tk,
    sequences: Mapping[str, str],
    *,
    color_bases: bool = True,
    show_grid: bool = False,
) -> None:
    """
    Fast virtualized FASTA/alignment viewer.

    Key optimization:
    - The canvas scrollregion represents the full alignment.
    - Only the currently visible rows/columns are rendered.
    - Sequence letters are drawn as one text string per visible row.
    - Optional base background colors are drawn only for visible cells.

    This is suitable for large alignments/sequences where the naive
    create-rectangle/create-text-per-base approach becomes extremely slow.
    """

    if not sequences:
        viewer = tk.Toplevel(root)
        viewer.title("FASTA Viewer")
        viewer.geometry("800x300")
        tk.Label(viewer, text="No sequences to display.").pack(expand=True)
        return

    headers = list(sequences.keys())
    seqs = [str(sequences[h]) for h in headers]

    n_rows = len(headers)
    seq_length = max(len(seq) for seq in seqs)

    viewer = tk.Toplevel(root)
    viewer.title("FASTA Viewer")
    viewer.geometry("1200x500")

    font = tkfont.Font(family="Courier New", size=10)
    bold_font = font.copy()
    bold_font.configure(weight="bold")

    name_font = tkfont.nametofont("TkDefaultFont")
    name_bold_font = name_font.copy()
    name_bold_font.configure(weight="bold")

    base_width = max(1, font.measure("W"))
    row_height = max(20, font.metrics("linespace") + 8)

    name_width = max(
        name_font.measure("Sequence name"),
        max(name_font.measure(header) for header in headers),
    ) + 16

    position_digits = len(str(seq_length))
    header_height = max(34, position_digits * 9 + 14)

    total_width = seq_length * base_width
    total_height = n_rows * row_height

    outer = tk.Frame(viewer)
    outer.pack(fill="both", expand=True)

    corner_canvas = tk.Canvas(
        outer,
        width=name_width,
        height=header_height,
        highlightthickness=0,
        background="white",
    )
    top_canvas = tk.Canvas(
        outer,
        height=header_height,
        highlightthickness=0,
        background="white",
    )
    name_canvas = tk.Canvas(
        outer,
        width=name_width,
        highlightthickness=0,
        background="white",
    )
    seq_canvas = tk.Canvas(
        outer,
        highlightthickness=0,
        background="white",
    )

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

    for canvas in (top_canvas, name_canvas, seq_canvas):
        canvas.configure(borderwidth=0)

    top_canvas.configure(scrollregion=(0, 0, total_width, header_height))
    name_canvas.configure(scrollregion=(0, 0, name_width, total_height))
    seq_canvas.configure(scrollregion=(0, 0, total_width, total_height))

    # Makes wheel/keyboard scrolling align naturally to residues/rows.
    seq_canvas.configure(xscrollincrement=base_width, yscrollincrement=row_height)
    top_canvas.configure(xscrollincrement=base_width)
    name_canvas.configure(yscrollincrement=row_height)

    hovered_row: int | None = None
    redraw_pending = False

    def request_redraw() -> None:
        nonlocal redraw_pending

        if redraw_pending:
            return

        redraw_pending = True

        def _do_redraw() -> None:
            nonlocal redraw_pending
            redraw_pending = False
            redraw()

        viewer.after_idle(_do_redraw)

    def sync_y_scrollbar(first: str, last: str) -> None:
        y_scroll.set(first, last)

    def sync_x_scrollbar(first: str, last: str) -> None:
        x_scroll.set(first, last)

    seq_canvas.configure(
        xscrollcommand=sync_x_scrollbar,
        yscrollcommand=sync_y_scrollbar,
    )

    def yview(*args: object) -> None:
        name_canvas.yview(*args)
        seq_canvas.yview(*args)
        request_redraw()

    def xview(*args: object) -> None:
        top_canvas.xview(*args)
        seq_canvas.xview(*args)
        request_redraw()

    y_scroll.configure(command=yview)
    x_scroll.configure(command=xview)

    def draw_corner() -> None:
        corner_canvas.delete("all")
        corner_canvas.create_rectangle(
            0,
            0,
            name_width,
            header_height,
            fill="white",
            outline="gray",
        )
        corner_canvas.create_text(
            8,
            header_height / 2,
            text="Sequence name",
            anchor="w",
            font=name_bold_font,
        )

    def visible_range() -> tuple[int, int, int, int]:
        x0 = max(0, int(seq_canvas.canvasx(0)))
        y0 = max(0, int(seq_canvas.canvasy(0)))

        width = max(1, seq_canvas.winfo_width())
        height = max(1, seq_canvas.winfo_height())

        first_col = max(0, x0 // base_width)
        last_col = min(seq_length, (x0 + width) // base_width + 2)

        first_row = max(0, y0 // row_height)
        last_row = min(n_rows, (y0 + height) // row_height + 2)

        return first_row, last_row, first_col, last_col

    def padded_slice(sequence: str, first_col: int, last_col: int) -> str:
        wanted = last_col - first_col

        if first_col >= len(sequence):
            return "-" * wanted

        segment = sequence[first_col:last_col]

        if len(segment) < wanted:
            segment += "-" * (wanted - len(segment))

        return segment

    def base_fill(base: str) -> str:
        return "#" + COLORS.get(base.upper(), "FFFFFF")

    def draw_header(first_col: int, last_col: int) -> None:
        top_canvas.delete("all")

        x_left = first_col * base_width
        x_right = last_col * base_width

        top_canvas.create_rectangle(
            x_left,
            0,
            x_right,
            header_height,
            fill="white",
            outline="",
        )

        # Draw a label every 10 columns, plus first visible column.
        # Drawing every single position label is possible but visually noisy.
        start = first_col + 1
        end = last_col

        first_labeled = start
        if first_labeled % 10 != 0:
            first_labeled += 10 - (first_labeled % 10)

        label_positions = {start}
        label_positions.update(range(first_labeled, end + 1, 10))

        for pos_1_based in sorted(label_positions):
            col = pos_1_based - 1
            x = col * base_width + base_width / 2

            top_canvas.create_line(
                x,
                header_height - 8,
                x,
                header_height,
                fill="gray",
            )
            top_canvas.create_text(
                x,
                header_height - 10,
                text=str(pos_1_based),
                anchor="s",
                font=name_font,
            )

        if show_grid:
            for col in range(first_col, last_col + 1):
                x = col * base_width
                top_canvas.create_line(x, 0, x, header_height, fill="#e6e6e6")

        top_canvas.create_line(
            x_left,
            header_height - 1,
            x_right,
            header_height - 1,
            fill="gray",
        )

    def draw_names(first_row: int, last_row: int) -> None:
        name_canvas.delete("all")

        for row_index in range(first_row, last_row):
            y1 = row_index * row_height
            y2 = y1 + row_height
            is_hovered = row_index == hovered_row

            fill = "#f2f2f2" if is_hovered else "white"
            font_to_use = name_bold_font if is_hovered else name_font

            name_canvas.create_rectangle(
                0,
                y1,
                name_width,
                y2,
                fill=fill,
                outline="#d0d0d0",
            )
            name_canvas.create_text(
                8,
                y1 + row_height / 2,
                text=headers[row_index],
                anchor="w",
                font=font_to_use,
            )

    def draw_sequences(
        first_row: int,
        last_row: int,
        first_col: int,
        last_col: int,
    ) -> None:
        seq_canvas.delete("all")

        x_start = first_col * base_width
        x_end = last_col * base_width

        for row_index in range(first_row, last_row):
            y1 = row_index * row_height
            y2 = y1 + row_height
            is_hovered = row_index == hovered_row

            visible_seq = padded_slice(seqs[row_index], first_col, last_col)

            if color_bases:
                for offset, base in enumerate(visible_seq):
                    x1 = (first_col + offset) * base_width
                    x2 = x1 + base_width

                    seq_canvas.create_rectangle(
                        x1,
                        y1,
                        x2,
                        y2,
                        fill=base_fill(base),
                        outline="#e6e6e6" if show_grid else "",
                    )
            else:
                seq_canvas.create_rectangle(
                    x_start,
                    y1,
                    x_end,
                    y2,
                    fill="#f2f2f2" if is_hovered else "white",
                    outline="",
                )

            if is_hovered:
                seq_canvas.create_rectangle(
                    x_start,
                    y1,
                    x_end,
                    y2,
                    outline="black",
                    width=1,
                )

            seq_canvas.create_text(
                x_start,
                y1 + 4,
                text=visible_seq,
                anchor="nw",
                font=bold_font if is_hovered else font,
                fill="black",
            )

        if show_grid:
            for row_index in range(first_row, last_row + 1):
                y = row_index * row_height
                seq_canvas.create_line(x_start, y, x_end, y, fill="#e6e6e6")

    def redraw() -> None:
        first_row, last_row, first_col, last_col = visible_range()

        draw_header(first_col, last_col)
        draw_names(first_row, last_row)
        draw_sequences(first_row, last_row, first_col, last_col)

    def row_from_event(canvas: tk.Canvas, event: tk.Event) -> int | None:
        y = int(canvas.canvasy(event.y))
        row_index = y // row_height

        if 0 <= row_index < n_rows:
            return row_index

        return None

    def on_row_motion(canvas: tk.Canvas, event: tk.Event) -> str:
        nonlocal hovered_row

        row_index = row_from_event(canvas, event)

        if row_index != hovered_row:
            hovered_row = row_index
            request_redraw()

        return "break"

    def clear_hover() -> None:
        nonlocal hovered_row

        if hovered_row is not None:
            hovered_row = None
            request_redraw()

    def scroll_vertical(units: int) -> None:
        name_canvas.yview_scroll(units, "units")
        seq_canvas.yview_scroll(units, "units")
        clear_hover()
        request_redraw()

    def scroll_horizontal(units: int) -> None:
        top_canvas.xview_scroll(units, "units")
        seq_canvas.xview_scroll(units, "units")
        clear_hover()
        request_redraw()

    def wheel_units(event: tk.Event) -> int:
        delta = getattr(event, "delta", 0)

        if delta == 0:
            return 0

        if abs(delta) >= 120:
            return -int(delta / 120)

        return -1 if delta > 0 else 1

    def wants_horizontal_scroll(event: tk.Event) -> bool:
        # Shift is standard horizontal scroll on many systems.
        # Ctrl is kept because your original implementation used it.
        shift_pressed = bool(event.state & 0x0001)
        ctrl_pressed = bool(event.state & 0x0004)
        return shift_pressed or ctrl_pressed

    def on_mousewheel(event: tk.Event) -> str:
        units = wheel_units(event)

        if units == 0:
            return "break"

        if wants_horizontal_scroll(event):
            scroll_horizontal(units * 3)
        else:
            scroll_vertical(units * 3)

        return "break"

    def on_linux_scroll_up(event: tk.Event) -> str:
        if wants_horizontal_scroll(event):
            scroll_horizontal(-3)
        else:
            scroll_vertical(-3)

        return "break"

    def on_linux_scroll_down(event: tk.Event) -> str:
        if wants_horizontal_scroll(event):
            scroll_horizontal(3)
        else:
            scroll_vertical(3)

        return "break"

    def on_configure(_event: tk.Event) -> None:
        request_redraw()

    name_canvas.bind("<Motion>", lambda event: on_row_motion(name_canvas, event))
    seq_canvas.bind("<Motion>", lambda event: on_row_motion(seq_canvas, event))

    for canvas in (corner_canvas, top_canvas, name_canvas, seq_canvas):
        canvas.bind("<MouseWheel>", on_mousewheel)
        canvas.bind("<Button-4>", on_linux_scroll_up)
        canvas.bind("<Button-5>", on_linux_scroll_down)
        canvas.bind("<Leave>", lambda _event: clear_hover())

    seq_canvas.bind("<Configure>", on_configure)
    name_canvas.bind("<Configure>", on_configure)
    top_canvas.bind("<Configure>", on_configure)

    viewer.bind("<Leave>", lambda _event: clear_hover())

    draw_corner()
    request_redraw()