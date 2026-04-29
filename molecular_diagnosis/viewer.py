import tkinter as tk
import tkinter.font as tkfont
from collections.abc import Mapping

from PIL import Image, ImageDraw, ImageTk

from molecular_diagnosis.constants import COLORS


def open_fasta_viewer(
    root: tk.Tk,
    sequences: Mapping[str, str],
    *,
    color_bases: bool = True,
    show_letters: bool = True,
    show_grid: bool = False,
) -> None:
    """
    Fast bitmap-backed FASTA/alignment viewer.

    This avoids creating Canvas rectangles/text for every visible residue.
    The sequence viewport is rendered into a single PIL image, then displayed
    as one Canvas image item.

    This is generally faster than a normal Tkinter Canvas-grid approach,
    especially during horizontal scrolling.
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

    # Tk font for measurement and Canvas text.
    seq_font = tkfont.Font(family="Courier New", size=10)
    name_font = tkfont.nametofont("TkDefaultFont")
    name_bold_font = name_font.copy()
    name_bold_font.configure(weight="bold")

    base_width = max(7, seq_font.measure("W"))
    row_height = max(20, seq_font.metrics("linespace") + 8)

    name_width = max(
        name_font.measure("Sequence name"),
        max(name_font.measure(header) for header in headers),
    ) + 16

    header_height = 42

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

    x_offset = 0
    y_offset = 0
    hovered_row: int | None = None
    redraw_pending = False

    selection_anchor: tuple[int, int] | None = None
    selection_active: tuple[int, int] | None = None

    # Keep a reference, otherwise Tk may garbage-collect the image.
    seq_canvas._viewport_image = None  # type: ignore[attr-defined]

    color_cache: dict[str, tuple[int, int, int]] = {}

    def parse_color(base: str) -> tuple[int, int, int]:
        base = base.upper()

        if base in color_cache:
            return color_cache[base]

        hex_color = COLORS.get(base, "FFFFFF").lstrip("#")

        if len(hex_color) != 6:
            rgb = (255, 255, 255)
        else:
            try:
                rgb = (
                    int(hex_color[0:2], 16),
                    int(hex_color[2:4], 16),
                    int(hex_color[4:6], 16),
                )
            except ValueError:
                rgb = (255, 255, 255)

        color_cache[base] = rgb
        return rgb

    def clamp_offsets() -> None:
        nonlocal x_offset, y_offset

        viewport_w = max(1, seq_canvas.winfo_width())
        viewport_h = max(1, seq_canvas.winfo_height())

        max_x = max(0, total_width - viewport_w)
        max_y = max(0, total_height - viewport_h)

        x_offset = min(max(0, x_offset), max_x)
        y_offset = min(max(0, y_offset), max_y)

    def update_scrollbars() -> None:
        viewport_w = max(1, seq_canvas.winfo_width())
        viewport_h = max(1, seq_canvas.winfo_height())

        if total_width <= 0:
            x_scroll.set(0, 1)
        else:
            first = x_offset / total_width
            last = min(1, (x_offset + viewport_w) / total_width)
            x_scroll.set(first, last)

        if total_height <= 0:
            y_scroll.set(0, 1)
        else:
            first = y_offset / total_height
            last = min(1, (y_offset + viewport_h) / total_height)
            y_scroll.set(first, last)

    def visible_range() -> tuple[int, int, int, int, int, int]:
        viewport_w = max(1, seq_canvas.winfo_width())
        viewport_h = max(1, seq_canvas.winfo_height())

        first_col = max(0, x_offset // base_width)
        last_col = min(seq_length, (x_offset + viewport_w) // base_width + 2)

        first_row = max(0, y_offset // row_height)
        last_row = min(n_rows, (y_offset + viewport_h) // row_height + 2)

        x_remainder = x_offset - first_col * base_width
        y_remainder = y_offset - first_row * row_height

        return first_row, last_row, first_col, last_col, x_remainder, y_remainder

    def padded_slice(sequence: str, first_col: int, last_col: int) -> str:
        wanted = last_col - first_col

        if first_col >= len(sequence):
            return "-" * wanted

        segment = sequence[first_col:last_col]

        if len(segment) < wanted:
            segment += "-" * (wanted - len(segment))

        return segment

    def request_redraw() -> None:
        nonlocal redraw_pending

        if redraw_pending:
            return

        redraw_pending = True

        def _redraw() -> None:
            nonlocal redraw_pending
            redraw_pending = False
            redraw()

        viewer.after_idle(_redraw)

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

    def draw_header(first_col: int, last_col: int, x_remainder: int) -> None:
        top_canvas.delete("all")

        viewport_w = max(1, top_canvas.winfo_width())

        top_canvas.create_rectangle(
            0,
            0,
            viewport_w,
            header_height,
            fill="white",
            outline="",
        )

        # Positions are 1-based for display.
        start_pos = first_col + 1
        end_pos = last_col

        # Minor tick every 5.
        first_minor = start_pos
        if first_minor % 5 != 0:
            first_minor += 5 - (first_minor % 5)

        for pos_1_based in range(first_minor, end_pos + 1, 5):
            col = pos_1_based - 1
            x = (col - first_col) * base_width - x_remainder + base_width / 2

            if 0 <= x <= viewport_w:
                top_canvas.create_line(
                    x,
                    header_height - 6,
                    x,
                    header_height,
                    fill="gray",
                )

        # Major numbered tick every 10.
        first_major = start_pos
        if first_major % 10 != 0:
            first_major += 10 - (first_major % 10)

        for pos_1_based in range(first_major, end_pos + 1, 10):
            col = pos_1_based - 1
            x = (col - first_col) * base_width - x_remainder + base_width / 2
            label = str(pos_1_based)

            label_width = name_font.measure(label)

            # Prevent clipped/half-visible numbers at the left and right edges.
            if label_width / 2 <= x <= viewport_w - label_width / 2:
                top_canvas.create_line(
                    x,
                    header_height - 12,
                    x,
                    header_height,
                    fill="black",
                )
                top_canvas.create_text(
                    x,
                    header_height - 14,
                    text=label,
                    anchor="s",
                    font=name_font,
                )

        top_canvas.create_line(
            0,
            header_height - 1,
            viewport_w,
            header_height - 1,
            fill="gray",
        )
        
    def draw_names(first_row: int, last_row: int, y_remainder: int) -> None:
        name_canvas.delete("all")

        viewport_h = max(1, name_canvas.winfo_height())

        name_canvas.create_rectangle(
            0,
            0,
            name_width,
            viewport_h,
            fill="white",
            outline="",
        )

        for row_index in range(first_row, last_row):
            y1 = (row_index - first_row) * row_height - y_remainder
            y2 = y1 + row_height

            is_hovered = row_index == hovered_row
            fill = "#f2f2f2" if is_hovered else "white"

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
                font=name_bold_font if is_hovered else name_font,
            )

    def selected_rectangle() -> tuple[int, int, int, int] | None:
        if selection_anchor is None or selection_active is None:
            return None

        row_a, col_a = selection_anchor
        row_b, col_b = selection_active

        row_start = min(row_a, row_b)
        row_end = max(row_a, row_b)
        col_start = min(col_a, col_b)
        col_end = max(col_a, col_b)

        return row_start, row_end, col_start, col_end

    def render_sequence_bitmap(
        first_row: int,
        last_row: int,
        first_col: int,
        last_col: int,
        x_remainder: int,
        y_remainder: int,
    ) -> ImageTk.PhotoImage:
        viewport_w = max(1, seq_canvas.winfo_width())
        viewport_h = max(1, seq_canvas.winfo_height())

        image = Image.new("RGB", (viewport_w, viewport_h), "white")
        draw = ImageDraw.Draw(image)

        for row_index in range(first_row, last_row):
            y1 = (row_index - first_row) * row_height - y_remainder
            y2 = y1 + row_height

            if y2 < 0 or y1 > viewport_h:
                continue

            is_hovered = row_index == hovered_row

            if not color_bases:
                row_fill = (242, 242, 242) if is_hovered else (255, 255, 255)
                draw.rectangle((0, y1, viewport_w, y2), fill=row_fill)
            else:
                visible_seq = padded_slice(seqs[row_index], first_col, last_col)

                for offset, base in enumerate(visible_seq):
                    x1 = offset * base_width - x_remainder
                    x2 = x1 + base_width

                    if x2 < 0 or x1 > viewport_w:
                        continue

                    draw.rectangle(
                        (x1, y1, x2, y2),
                        fill=parse_color(base),
                    )

                    if show_grid:
                        draw.rectangle(
                            (x1, y1, x2, y2),
                            outline=(230, 230, 230),
                        )

            if is_hovered:
                draw.rectangle(
                    (0, y1, viewport_w, y2),
                    outline=(0, 0, 0),
                )
            
            selected = selected_rectangle()

            if selected is not None:
                row_start, row_end, col_start, col_end = selected

                if row_start <= row_index <= row_end:
                    visible_col_start = max(first_col, col_start)
                    visible_col_end = min(last_col - 1, col_end)

                    if visible_col_start <= visible_col_end:
                        sx1 = (
                            (visible_col_start - first_col) * base_width
                            - x_remainder
                        )
                        sx2 = (
                            (visible_col_end - first_col + 1) * base_width
                            - x_remainder
                        )

                        draw.rectangle(
                            (sx1, y1, sx2, y2),
                            outline=(0, 0, 0),
                            width=2,
                        )

        return ImageTk.PhotoImage(image)

    def draw_sequence_letters(
        first_row: int,
        last_row: int,
        first_col: int,
        last_col: int,
        x_remainder: int,
        y_remainder: int,
    ) -> None:
        if not show_letters:
            return

        for row_index in range(first_row, last_row):
            y = (row_index - first_row) * row_height - y_remainder + 4
            x = -x_remainder

            visible_seq = padded_slice(seqs[row_index], first_col, last_col)

            seq_canvas.create_text(
                x,
                y,
                text=visible_seq,
                anchor="nw",
                font=seq_font,
                fill="black",
            )

    def redraw() -> None:
        clamp_offsets()
        update_scrollbars()

        (
            first_row,
            last_row,
            first_col,
            last_col,
            x_remainder,
            y_remainder,
        ) = visible_range()

        draw_header(first_col, last_col, x_remainder)
        draw_names(first_row, last_row, y_remainder)

        seq_canvas.delete("all")

        photo = render_sequence_bitmap(
            first_row,
            last_row,
            first_col,
            last_col,
            x_remainder,
            y_remainder,
        )

        seq_canvas._viewport_image = photo  # type: ignore[attr-defined]
        seq_canvas.create_image(0, 0, image=photo, anchor="nw")

        # Text is still Canvas text, but only one item per visible row,
        # not one item per base.
        draw_sequence_letters(
            first_row,
            last_row,
            first_col,
            last_col,
            x_remainder,
            y_remainder,
        )

    def xview(*args: object) -> None:
        nonlocal x_offset

        if not args:
            return

        command = args[0]

        if command == "moveto":
            fraction = float(args[1])
            x_offset = int(fraction * total_width)

        elif command == "scroll":
            amount = int(args[1])
            unit = str(args[2])

            if unit == "pages":
                x_offset += amount * max(1, seq_canvas.winfo_width() - base_width)
            else:
                x_offset += amount * base_width * 5

        request_redraw()

    def yview(*args: object) -> None:
        nonlocal y_offset

        if not args:
            return

        command = args[0]

        if command == "moveto":
            fraction = float(args[1])
            y_offset = int(fraction * total_height)

        elif command == "scroll":
            amount = int(args[1])
            unit = str(args[2])

            if unit == "pages":
                y_offset += amount * max(1, seq_canvas.winfo_height() - row_height)
            else:
                y_offset += amount * row_height * 3

        request_redraw()

    x_scroll.configure(command=xview)
    y_scroll.configure(command=yview)

    def row_from_event(event: tk.Event) -> int | None:
        row_index = (y_offset + event.y) // row_height

        if 0 <= row_index < n_rows:
            return row_index

        return None

    def on_motion(event: tk.Event) -> str:
        nonlocal hovered_row

        row_index = row_from_event(event)

        if row_index != hovered_row:
            hovered_row = row_index
            request_redraw()

        return "break"

    def clear_hover() -> None:
        nonlocal hovered_row

        if hovered_row is not None:
            hovered_row = None
            request_redraw()

    def wants_horizontal_scroll(event: tk.Event) -> bool:
        shift_pressed = bool(event.state & 0x0001)
        ctrl_pressed = bool(event.state & 0x0004)
        return shift_pressed or ctrl_pressed

    def wheel_units(event: tk.Event) -> int:
        delta = getattr(event, "delta", 0)

        if delta == 0:
            return 0

        if abs(delta) >= 120:
            return -int(delta / 120)

        return -1 if delta > 0 else 1

    def on_mousewheel(event: tk.Event) -> str:
        nonlocal x_offset, y_offset

        units = wheel_units(event)

        if units == 0:
            return "break"

        if wants_horizontal_scroll(event):
            x_offset += units * base_width * 8
        else:
            y_offset += units * row_height * 3

        request_redraw()
        return "break"

    def on_linux_scroll_up(event: tk.Event) -> str:
        nonlocal x_offset, y_offset

        if wants_horizontal_scroll(event):
            x_offset -= base_width * 8
        else:
            y_offset -= row_height * 3

        request_redraw()
        return "break"

    def on_linux_scroll_down(event: tk.Event) -> str:
        nonlocal x_offset, y_offset

        if wants_horizontal_scroll(event):
            x_offset += base_width * 8
        else:
            y_offset += row_height * 3

        request_redraw()
        return "break"

    def on_configure(_event: tk.Event) -> None:
        request_redraw()

    seq_canvas.bind("<Motion>", on_motion)
    name_canvas.bind("<Motion>", on_motion)

    for canvas in (corner_canvas, top_canvas, name_canvas, seq_canvas):
        canvas.bind("<MouseWheel>", on_mousewheel)
        canvas.bind("<Button-4>", on_linux_scroll_up)
        canvas.bind("<Button-5>", on_linux_scroll_down)
        canvas.bind("<Leave>", lambda _event: clear_hover())

    seq_canvas.bind("<Configure>", on_configure)
    name_canvas.bind("<Configure>", on_configure)
    top_canvas.bind("<Configure>", on_configure)

    viewer.bind("<Leave>", lambda _event: clear_hover())

    def cell_from_event(event: tk.Event) -> tuple[int, int] | None:
        row = (y_offset + event.y) // row_height
        col = (x_offset + event.x) // base_width

        if 0 <= row < n_rows and 0 <= col < seq_length:
            return int(row), int(col)

        return None

    def on_press(event: tk.Event) -> str:
        nonlocal selection_anchor, selection_active

        cell = cell_from_event(event)

        if cell is not None:
            selection_anchor = cell
            selection_active = cell

            row, col = cell
            header = headers[row]
            base = seqs[row][col] if col < len(seqs[row]) else "-"

            print(f"Clicked row={row}, header={header}, col={col + 1}, base={base}")

            request_redraw()

        return "break"

    def on_drag(event: tk.Event) -> str:
        nonlocal selection_active

        cell = cell_from_event(event)

        if cell is not None and cell != selection_active:
            selection_active = cell
            request_redraw()

        return "break"

    def on_release(event: tk.Event) -> str:
        return "break"

    seq_canvas.bind("<Button-1>", on_press)
    seq_canvas.bind("<B1-Motion>", on_drag)
    seq_canvas.bind("<ButtonRelease-1>", on_release)

    draw_corner()
    request_redraw() 