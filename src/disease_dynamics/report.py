"""Status lines for notebooks: PASS / WARN / FAIL / INFO, collected into one panel.

The log only records; the notebook decides how to display it. ``raise_if_failed`` is called
in the last cell, so a notebook with a failed check cannot finish quietly.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class StatusLog:
    """Ordered record of checks, grouped by section."""

    entries: list[tuple[str, str, str]] = field(default_factory=list)
    section: str = ""

    def start(self, section: str) -> None:
        """Set the section name used for the following entries."""
        self.section = section

    def check(self, condition: bool, message: str) -> str:
        """Record PASS or FAIL and return the formatted line."""
        return self._add("PASS" if condition else "FAIL", message)

    def warn(self, message: str) -> str:
        """Record a warning: not wrong, but the reader must know."""
        return self._add("WARN", message)

    def info(self, message: str) -> str:
        """Record a result that is reported, not asserted."""
        return self._add("INFO", message)

    def _add(self, level: str, message: str) -> str:
        self.entries.append((self.section, level, message))
        return f"{level:4}  {message}"

    def counts(self) -> dict[str, int]:
        """Number of entries per level."""
        result = {"PASS": 0, "WARN": 0, "FAIL": 0, "INFO": 0}
        for _, level, _ in self.entries:
            result[level] += 1
        return result

    def markdown(self) -> str:
        """All entries as a Markdown table."""
        rows = ["| Section | Status | Check |", "|---|---|---|"]
        rows += [f"| {s} | {lvl} | {m} |" for s, lvl, m in self.entries]
        return "\n".join(rows)

    def raise_if_failed(self) -> None:
        """Raise AssertionError listing every failed check."""
        failed = [f"[{s}] {m}" for s, lvl, m in self.entries if lvl == "FAIL"]
        if failed:
            raise AssertionError(f"{len(failed)} check(s) failed:\n" + "\n".join(failed))
