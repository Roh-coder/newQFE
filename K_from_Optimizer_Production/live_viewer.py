"""
live_viewer.py — minimal live PNG viewer for the optimizer frame.

Polls a single rolling PNG (written by ``visualization.OptimizerPlotter``)
and re-displays it in a Tk window every ``poll_ms`` milliseconds.

Designed to coexist with the optimizer running in a background thread:
the Tk mainloop owns the *main* thread; the optimizer is started via
``start_optimizer_thread(target_callable)`` and the window closes
automatically when the optimizer finishes.

Stdlib only — no extra deps beyond what matplotlib already pulls in.
Tk's PhotoImage natively reads PNG on Tk ≥ 8.6 (Python ≥ 3.0 ships 8.6
on every supported OS).
"""
from __future__ import annotations

import os
import threading
import time
from typing import Callable, Optional


class LiveViewer:
    """Tk window that reloads ``png_path`` whenever its mtime changes."""

    def __init__(self, png_path: str, *,
                 title: str = "K_from_Optimizer — live",
                 poll_ms: int = 500,
                 max_width: int = 1400,
                 max_height: int = 900):
        self.png_path  = png_path
        self.poll_ms   = int(poll_ms)
        self.max_w     = int(max_width)
        self.max_h     = int(max_height)
        self._last_mtime: Optional[float] = None
        self._optimizer_thread: Optional[threading.Thread] = None
        self._optimizer_done = threading.Event()
        self._exception: Optional[BaseException] = None

        import tkinter as tk
        self._tk = tk
        self.root = tk.Tk()
        self.root.title(title)
        self.label = tk.Label(self.root, text="Waiting for first frame…",
                              compound="center")
        self.label.pack(fill="both", expand=True)
        self.status = tk.Label(self.root,
                               text=f"watching {os.path.basename(png_path)}",
                               anchor="w")
        self.status.pack(fill="x")
        # Hold the most recent PhotoImage so Tk doesn't garbage-collect it.
        self._photo = None
        self.root.protocol("WM_DELETE_WINDOW", self._on_close)

    # ------------------------------------------------------------------
    def start_optimizer_thread(self, target: Callable[[], object]) -> None:
        """Launch ``target()`` in a daemon thread; viewer closes when done."""
        def _runner():
            try:
                target()
            except BaseException as e:   # surface to main thread on close
                self._exception = e
            finally:
                self._optimizer_done.set()

        self._optimizer_thread = threading.Thread(
            target=_runner, name="optimizer", daemon=True)
        self._optimizer_thread.start()

    # ------------------------------------------------------------------
    def _on_close(self) -> None:
        # Tk window closed by the user — the optimizer thread (daemon)
        # will be killed when the process exits.  We just stop the
        # mainloop here.
        try:
            self.root.destroy()
        except Exception:
            pass

    def _refresh(self) -> None:
        # Stop polling once the optimizer has finished AND the user has had
        # a chance to see the final frame (we keep the window open until
        # they close it).
        try:
            if os.path.exists(self.png_path):
                m = os.path.getmtime(self.png_path)
                if self._last_mtime is None or m > self._last_mtime:
                    self._last_mtime = m
                    self._reload_image()
            done_msg = ""
            if self._optimizer_done.is_set():
                if self._exception is None:
                    done_msg = "  [optimizer DONE — close window to exit]"
                else:
                    done_msg = f"  [optimizer FAILED: {self._exception}]"
            self.status.configure(
                text=f"{os.path.basename(self.png_path)}"
                     f"   mtime={time.strftime('%H:%M:%S', time.localtime(self._last_mtime)) if self._last_mtime else '—'}"
                     f"{done_msg}")
        finally:
            self.root.after(self.poll_ms, self._refresh)

    def _reload_image(self) -> None:
        # Use Pillow if present (handles arbitrary PNG depths) but fall
        # back to Tk's PhotoImage which is stdlib-only and works for
        # standard 8-bit/channel matplotlib output.
        try:
            from PIL import Image, ImageTk
            img = Image.open(self.png_path)
            img.thumbnail((self.max_w, self.max_h))
            self._photo = ImageTk.PhotoImage(img)
        except Exception:
            try:
                self._photo = self._tk.PhotoImage(file=self.png_path)
            except Exception as e:
                self.label.configure(text=f"<failed to load PNG: {e}>",
                                     image="")
                return
        self.label.configure(image=self._photo, text="")

    # ------------------------------------------------------------------
    def mainloop(self) -> None:
        self.root.after(self.poll_ms, self._refresh)
        self.root.mainloop()
        # If the optimizer thread raised, re-raise on the main thread so
        # the process exits with a non-zero code.
        if self._exception is not None:
            raise self._exception
