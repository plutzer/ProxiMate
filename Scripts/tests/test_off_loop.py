"""off_loop runs a slow computation in a worker thread and delivers only the latest run."""

import asyncio
import threading

from shiny import reactive

from off_loop import off_loop


async def _watch(compute, until):
    """Record each value ``compute()`` delivers to an effect, flushing until ``until(values)``."""
    delivered = []

    @reactive.effect
    def record():
        delivered.append(compute())

    for _ in range(300):
        await reactive.flush()
        if until(delivered):
            return delivered, record
        await asyncio.sleep(0.01)
    raise AssertionError(f"no expected result; delivered {delivered}")


def test_runs_in_a_worker_thread():
    def work(x):
        return x * 10, threading.get_ident()

    async def run():
        value = reactive.value(4)
        delivered, _ = await _watch(off_loop(work, lambda: (value(),)), lambda d: d)
        return delivered

    [(result, worker)] = asyncio.run(run())
    assert result == 40
    assert worker != threading.get_ident()


def test_delivering_a_result_does_not_start_another_run():
    calls = []

    def work(x):
        calls.append(x)
        return x

    async def run():
        value = reactive.value(1)
        await _watch(off_loop(work, lambda: (value(),)), lambda d: d)
        for _ in range(10):
            await asyncio.sleep(0.01)
            await reactive.flush()

    asyncio.run(run())
    assert calls == [1]


def test_new_arguments_replace_a_run_in_progress():
    release = threading.Event()
    calls = []

    def work(x):
        calls.append(x)
        if x == 1:
            release.wait(5)
        return x

    async def run():
        value = reactive.value(1)
        compute = off_loop(work, lambda: (value(),))
        delivered = []

        @reactive.effect
        def record():
            delivered.append(compute())

        await reactive.flush()
        await asyncio.sleep(0.05)
        value.set(2)
        await reactive.flush()
        release.set()
        for _ in range(30):
            await asyncio.sleep(0.01)
            await reactive.flush()
        return delivered

    assert asyncio.run(run()) == [2]
    assert calls == [1, 2]
