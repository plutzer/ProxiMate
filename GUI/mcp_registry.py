"""The operation registry behind the MCP ``call_tool`` tool.

An agent sees four MCP tools; everything it can do is an *operation* registered here
with a name, a mode, a one-line summary, tags for search, and a JSON schema derived
from the function's signature.  ``call`` validates arguments against that schema
before the function runs, so a misspelled or missing argument fails at once with a
message naming it, and every call is appended to an activity log with its actor.

Modes:

- ``read``      reads state; changes nothing
- ``sandbox``   computes from a dataset directory and returns the result; writes
                nothing under the dataset and changes nothing the GUI shows
- ``dataset``   creates or changes a dataset; the GUI's dataset table updates
- ``cytoscape`` acts on the network drawn in Cytoscape; the Cytoscape tab updates
"""

import collections
import dataclasses
import datetime
import inspect
import threading
import typing

MODES = {
    'read': "Reads state; changes nothing.",
    'sandbox': "Computes and returns a result. Writes nothing under the dataset and "
               "changes nothing the GUI shows; the GUI's own settings stay as the user left them.",
    'dataset': "Creates or changes a dataset directory. The GUI's dataset table and "
               "dropdowns update within a second; the change is recorded in the dataset's run.json.",
    'cytoscape': "Acts on the network drawn in the user's Cytoscape window. The Cytoscape tab's "
                 "status and activity panel update within a second, labeled with actor 'mcp'.",
}

_TYPES = {str: 'string', int: 'integer', float: 'number', bool: 'boolean',
          dict: 'object', list: 'array'}


@dataclasses.dataclass
class Op:
    name: str
    fn: typing.Callable
    mode: str
    summary: str
    tags: tuple
    schema: dict

    @property
    def doc(self):
        return inspect.getdoc(self.fn) or self.summary


OPS = {}
_ACTIVITY = collections.deque(maxlen=300)
_ACTIVITY_LOCK = threading.Lock()
_SEQ = 0          # number of calls logged so far; each entry carries its own value


def _json_type(annotation):
    """The JSON schema fragment for a parameter annotation.

    Supports the plain types, ``list[...]``/``dict[...]`` and ``Optional``; anything
    else is a registration error, so an op cannot silently take an unvalidated type.
    """
    origin = typing.get_origin(annotation)
    if origin is typing.Union or (origin is not None and origin.__name__ == 'UnionType'):
        members = [a for a in typing.get_args(annotation) if a is not type(None)]
        if len(members) != 1:
            raise TypeError(f"unsupported union annotation {annotation!r}")
        return {**_json_type(members[0]), 'nullable': True}
    if origin in (list, dict):
        return {'type': _TYPES[origin]}
    if annotation in _TYPES:
        return {'type': _TYPES[annotation]}
    raise TypeError(f"unsupported parameter annotation {annotation!r}")


def _schema(fn):
    sig = inspect.signature(fn)
    hints = typing.get_type_hints(fn)
    properties, required = {}, []
    for name, param in sig.parameters.items():
        if name not in hints:
            raise TypeError(f"{fn.__name__}: parameter {name!r} needs a type annotation")
        prop = _json_type(hints[name])
        if param.default is inspect.Parameter.empty:
            required.append(name)
        else:
            prop['default'] = param.default
        properties[name] = prop
    return {'type': 'object', 'properties': properties, 'required': required}


def register(mode, summary, tags=()):
    """Decorator: add a function to the registry under its own name."""
    if mode not in MODES:
        raise ValueError(f"mode must be one of {sorted(MODES)}")

    def wrap(fn):
        OPS[fn.__name__] = Op(fn.__name__, fn, mode, summary, tuple(tags), _schema(fn))
        return fn
    return wrap


def _get(name):
    if name not in OPS:
        raise KeyError(f"no operation named {name!r}; search_tools lists what exists")
    return OPS[name]


def search(query='', mode=None):
    """Operations whose name, summary or tags contain ``query`` (case-insensitive)."""
    if mode is not None and mode not in MODES:
        raise ValueError(f"mode must be one of {sorted(MODES)}")
    q = (query or '').strip().lower()
    hits = []
    for op in OPS.values():
        if mode and op.mode != mode:
            continue
        haystack = ' '.join([op.name, op.summary, *op.tags]).lower()
        if q and q not in haystack:
            continue
        hits.append({'name': op.name, 'mode': op.mode, 'summary': op.summary})
    return hits


def details(name):
    op = _get(name)
    return {'name': op.name, 'mode': op.mode, 'summary': op.summary, 'tags': list(op.tags),
            'effects': MODES[op.mode], 'doc': op.doc, 'schema': op.schema}


_CHECKS = {'string': lambda v: isinstance(v, str),
           'integer': lambda v: isinstance(v, int) and not isinstance(v, bool),
           'number': lambda v: isinstance(v, (int, float)) and not isinstance(v, bool),
           'boolean': lambda v: isinstance(v, bool),
           'object': lambda v: isinstance(v, dict),
           'array': lambda v: isinstance(v, list)}


def validate(name, arguments):
    """The arguments ``name`` accepts, checked against its schema; raises ValueError."""
    op = _get(name)
    arguments = dict(arguments or {})
    props = op.schema['properties']
    unexpected = [k for k in arguments if k not in props]
    if unexpected:
        raise ValueError(f"{name}: unexpected argument(s) {unexpected}; it takes {list(props)}")
    missing = [k for k in op.schema['required'] if k not in arguments]
    if missing:
        raise ValueError(f"{name}: missing required argument(s) {missing}")
    for key, value in arguments.items():
        prop = props[key]
        if value is None and prop.get('nullable'):
            continue
        if not _CHECKS[prop['type']](value):
            raise ValueError(f"{name}: argument {key!r} must be {prop['type']}, got {value!r}")
    return arguments


def call(name, arguments, actor='mcp'):
    """Validate, run and log one operation; exceptions propagate to the caller."""
    global _SEQ
    arguments = validate(name, arguments)
    op = _get(name)
    target = str(arguments.get('name') or arguments.get('dataset') or '')
    entry = {'ts': datetime.datetime.now().isoformat(timespec='seconds'), 'actor': actor,
             'op': name, 'mode': op.mode, 'ok': None, 'detail': target}
    try:
        result = op.fn(**arguments)
    except Exception as e:
        entry.update(ok=False, detail=f"{target}: {type(e).__name__}: {e}"[:300].strip(': '))
        raise
    else:
        entry['ok'] = True
        return result
    finally:
        with _ACTIVITY_LOCK:
            _SEQ += 1
            entry['seq'] = _SEQ
            _ACTIVITY.append(entry)


def activity(n_last=50):
    with _ACTIVITY_LOCK:
        return list(_ACTIVITY)[-n_last:]


def last_seq():
    return _SEQ


def activity_since(seq):
    """Entries logged after the call numbered ``seq``, oldest first."""
    with _ACTIVITY_LOCK:
        return [e for e in _ACTIVITY if e['seq'] > seq]
