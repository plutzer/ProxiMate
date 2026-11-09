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
    running_args = None

    @reactive.extended_task
    async def task(*fn_args):
        return await asyncio.to_thread(fn, *fn_args)

    def compute():
        nonlocal running_args
        fn_args = args()
        # args() hands back the same tuple object until something it reads changes, so a
        # different object means new arguments: stop the old run and start a new one.
        if fn_args is not running_args:
            running_args = fn_args
            task.cancel()
            task.invoke(*fn_args)
        return task.result()

    return compute
