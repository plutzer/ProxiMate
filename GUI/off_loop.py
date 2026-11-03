"""Slow computations for the GUI that run in a worker thread, off the event loop."""

import asyncio

from shiny import reactive


def off_loop(fn, args):
    """``compute()``, to call in a render or calc, returning ``fn(*args())`` computed in
    a worker thread.

    The event loop, and so every session and the websocket keepalive, stays free while
    ``fn`` runs; an output calling ``compute()`` shows as recalculating until it
    finishes.  ``args`` is a reactive calc, so ``fn`` reruns only when what ``args``
    reads changes, and only when something that calls ``compute()`` is shown.  New
    arguments drop the run in progress; only the latest arguments' result is delivered.
    """
    args = reactive.calc(args)
    started = [None]

    @reactive.extended_task
    async def task(*a):
        return await asyncio.to_thread(fn, *a)

    def compute():
        a = args()
        # The calc returns the same tuple until its inputs change; a new one is a new run.
        if a is not started[0]:
            started[0] = a
            task.cancel()
            task.invoke(*a)
        return task.result()

    return compute
