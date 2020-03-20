# Aubio Rust port via c2rust

Porting some Aubio library modules in pure (unsafe) Rust via __c2rust__.
The goal is to be as standalone as it can be (no libc / no static or dymamic lib)

The transpiled code is ugly unsafe Rust but it should work like the origin lib.
There is more idiomatic rust wrapper around.

### replace libc types
sed -i -- 's/libc::c_ulong/u64/g' ./src/transpiled/biquad.rs
sed -i -- 's/libc::c_float/f32/g' ./src/transpiled/biquad.rs
sed -i -- 's/libc::c_double/f64/g' ./src/transpiled/biquad.rs
sed -i -- 's/libc::c_void/core::ffi::c_void/g' ./src/transpiled/biquad.rs
sed -i -- 's/libc::c_uint/u32/g' ./src/transpiled/biquad.rs
sed -i -- 's/libc::c_int/i32/g' ./src/transpiled/biquad.rs
sed -i -- 's/libc::c_char/i8/g' ./src/transpiled/biquad.rs
sed -i -- 's/libc::c_uchar/u8/g' ./src/transpiled/biquad.rs
sed -i -- 's/libc::c_short/i16/g' ./src/transpiled/biquad.rs
sed -i -- 's/libc::c_longlong/i64/g' ./src/transpiled/biquad.rs