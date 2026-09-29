"""Copy the single shared workflow into an isolated work directory."""

from pathlib import Path

WORKFLOW = Path(__file__).resolve().parent / "workflow"


def install_workflow(contrast_id: str, work: Path) -> Path:
    """Install shared code once per comparison without editing the source tree."""
    code = work / "code" / contrast_id
    code.mkdir(parents=True, exist_ok=True)
    for source_file in WORKFLOW.rglob("*"):
        if not source_file.is_file() or "__pycache__" in source_file.parts:
            continue
        destination = code / source_file.relative_to(WORKFLOW)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(source_file.read_text(encoding="utf8"), encoding="utf8")
    return code
