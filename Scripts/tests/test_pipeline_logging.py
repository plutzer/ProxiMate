"""Static checks that pipeline modules report through the logger.

Output written with ``print`` goes to stdout, which is not carried by the
handlers that write ``proximate.log`` and, for a stage the GUI runs as a
subprocess, is not surfaced at the default log level either.  These checks are
static so they cost nothing and do not need fitted data to run.
"""

import ast
import os

import pytest

SCRIPTS_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO_DIR = os.path.dirname(SCRIPTS_DIR)
APP_PY = os.path.join(REPO_DIR, "GUI", "app.py")

# Modules that run as part of a scoring pipeline, where stray stdout is lost.
PIPELINE_MODULES = [
    "parse.py",
    "score.py",
    "annotator.py",
    "refactored_aft.py",
    "one_component_aft.py",
    "aft_impute_saint.py",
    "protein_groups.py",
    "compPASS_pval.py",
    "experimental_design.py",
    "bfdr_pool.py",
]


def _parse(filename):
    path = os.path.join(SCRIPTS_DIR, filename)
    with open(path, encoding="utf-8") as handle:
        return ast.parse(handle.read(), filename=path)


def _print_call_lines(tree):
    return [node.lineno for node in ast.walk(tree)
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id == "print"]


@pytest.mark.parametrize("filename", PIPELINE_MODULES)
def test_pipeline_modules_do_not_print(filename):
    lines = _print_call_lines(_parse(filename))

    assert not lines, f"{filename} calls print() at line(s) {lines}; use the logger"


@pytest.mark.parametrize("filename", PIPELINE_MODULES)
def test_pipeline_modules_use_the_shared_logger(filename):
    """Every pipeline module obtains its logger from log_config, so its records
    reach the dataset log and carry the run ID."""
    tree = _parse(filename)
    imports_log_config = any(
        isinstance(node, ast.ImportFrom) and node.module == "log_config"
        for node in ast.walk(tree))

    assert imports_log_config, f"{filename} does not import from log_config"


# ---------------------------------------------------------------------------
# GUI invariants
#
# GUI/app.py cannot be imported here (it needs shiny, which the analysis
# environment does not carry), so these read it as source instead.
# ---------------------------------------------------------------------------

def _app_tree():
    with open(APP_PY, encoding="utf-8") as handle:
        return ast.parse(handle.read(), filename=APP_PY)


def _calls_to(tree, attribute):
    return [node for node in ast.walk(tree)
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr == attribute]


def test_notifications_go_through_the_notify_helper():
    """A notification shown to the user must also be logged.  `notify` does
    both, so it has to be the only caller of ui.notification_show."""
    tree = _app_tree()
    notify_def = next(node for node in ast.walk(tree)
                      if isinstance(node, ast.FunctionDef) and node.name == "notify")
    inside_notify = {id(call) for call in _calls_to(notify_def, "notification_show")}

    stray = [call.lineno for call in _calls_to(tree, "notification_show")
             if id(call) not in inside_notify]

    assert not stray, (
        f"GUI/app.py calls ui.notification_show directly at line(s) {stray}; "
        "use notify() so the message is logged too")


def test_app_has_no_bare_except():
    """A bare `except:` also swallows KeyboardInterrupt and SystemExit."""
    tree = _app_tree()

    bare = [node.lineno for node in ast.walk(tree)
            if isinstance(node, ast.ExceptHandler) and node.type is None]

    assert not bare, f"GUI/app.py has a bare except at line(s) {bare}"


def test_app_does_not_print_tracebacks_to_the_console():
    """An end user of the containerized app cannot read the console; a traceback
    has to go to the logger to be recoverable."""
    tree = _app_tree()

    printed = [call.lineno for call in _calls_to(tree, "print_exc")]

    assert not printed, (
        f"GUI/app.py calls traceback.print_exc at line(s) {printed}; "
        "use logger.exception")
