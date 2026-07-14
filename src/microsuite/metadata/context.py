"""Workflow identity for stage-result envelopes.

Read once per ``stage_execution`` (never a process-global cache), so one process
running several datasets/benchmarks stays correct. An explicit ``WorkflowContext``
passed to ``stage_execution`` overrides the environment fallback.
"""

from __future__ import annotations

import os
from dataclasses import dataclass


@dataclass(frozen=True)
class WorkflowContext:
    run_id: str | None = None
    workflow_id: str | None = None
    workflow_run_id: str | None = None
    dataset_id: str | None = None


def workflow_context_from_env() -> WorkflowContext:
    """Build a WorkflowContext from ``MICROSUITE_*`` environment variables now."""
    return WorkflowContext(
        run_id=os.environ.get("MICROSUITE_RUN_ID"),
        workflow_id=os.environ.get("MICROSUITE_WORKFLOW_ID"),
        workflow_run_id=os.environ.get("MICROSUITE_WORKFLOW_RUN_ID"),
        dataset_id=os.environ.get("MICROSUITE_DATASET_ID"),
    )
