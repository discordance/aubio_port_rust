# Aubio Rust port via c2rust

Porting some Aubio library modules in pure (unsafe) Rust via __c2rust__.
The goal is to be as standalone as it can be (no libc / no static or dymamic lib)

The transpiled code is ugly unsafe Rust but it should work like the origin lib.
We try to provide an idiomatic rust wrapper around the libs.