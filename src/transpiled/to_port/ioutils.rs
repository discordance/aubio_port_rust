use ::libc;
extern "C" {
    #[no_mangle]
    fn aubio_log(level: sint_t, fmt: *const char_t, _: ...) -> uint_t;
    #[no_mangle]
    fn memcpy(_: *mut libc::c_void, _: *const libc::c_void, _: libc::c_ulong)
     -> *mut libc::c_void;
    #[no_mangle]
    fn memset(_: *mut libc::c_void, _: libc::c_int, _: libc::c_ulong)
     -> *mut libc::c_void;
}
/*
  Copyright (C) 2003-2015 Paul Brossier <piem@aubio.org>

  This file is part of aubio.

  aubio is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  aubio is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with aubio.  If not, see <http://www.gnu.org/licenses/>.

*/
/* * \file

  Definition of data types used in aubio

*/
/* * defined to 1 if aubio is compiled in double precision */
/* * short sample format (32 or 64 bits) */
pub type smpl_t = libc::c_float;
/* * print format for sample in double precision */
/* * unsigned integer */
pub type uint_t = libc::c_uint;
/* * signed integer */
pub type sint_t = libc::c_int;
/* * character */
pub type char_t = libc::c_char;
/*
  Copyright (C) 2003-2015 Paul Brossier <piem@aubio.org>

  This file is part of aubio.

  aubio is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  aubio is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with aubio.  If not, see <http://www.gnu.org/licenses/>.

*/
/* * \file

  Vector of real-valued data

  This file specifies the ::fvec_t buffer type, which is used throughout aubio
  to store vector of real-valued ::smpl_t.

  \example test-fvec.c

*/
/* * Buffer for real data

  Vector of real-valued data

  ::fvec_t is is the structure used to store vector of real-valued data, ::smpl_t .

  \code

  uint_t buffer_size = 1024;

  // create a vector of 512 values
  fvec_t * input = new_fvec (buffer_size);

  // set some values of the vector
  input->data[23] = 2.;
  // ..

  // compute the mean of the vector
  mean = fvec_mean(a_vector);

  // destroy the vector
  del_fvec(a_vector);

  \endcode

  See `examples/` and `tests/src` directories for more examples.

 */
#[derive(Copy, Clone)]
#[repr(C)]
pub struct fvec_t {
    pub length: uint_t,
    pub data: *mut smpl_t,
}
pub type C2RustUnnamed = libc::c_uint;
pub const AUBIO_FAIL: C2RustUnnamed = 1;
pub const AUBIO_OK: C2RustUnnamed = 0;
pub type aubio_log_level = libc::c_uint;
pub const AUBIO_LOG_LAST_LEVEL: aubio_log_level = 5;
pub const AUBIO_LOG_WRN: aubio_log_level = 4;
pub const AUBIO_LOG_DBG: aubio_log_level = 3;
pub const AUBIO_LOG_MSG: aubio_log_level = 2;
pub const AUBIO_LOG_INF: aubio_log_level = 1;
pub const AUBIO_LOG_ERR: aubio_log_level = 0;
#[derive(Copy, Clone)]
#[repr(C)]
pub struct fmat_t {
    pub length: uint_t,
    pub height: uint_t,
    pub data: *mut *mut smpl_t,
}
/*
  Copyright (C) 2016 Paul Brossier <piem@aubio.org>

  This file is part of aubio.

  aubio is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  aubio is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with aubio.  If not, see <http://www.gnu.org/licenses/>.

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_io_validate_samplerate(mut kind: *const char_t,
                                                      mut path: *const char_t,
                                                      mut samplerate: uint_t)
 -> uint_t {
    if samplerate as sint_t <= 0 as libc::c_int {
        aubio_log(AUBIO_LOG_ERR as libc::c_int,
                  b"AUBIO ERROR: %s: failed creating %s, samplerate should be positive, not %d\n\x00"
                      as *const u8 as *const libc::c_char, kind, path,
                  samplerate);
        return AUBIO_FAIL as libc::c_int as uint_t
    }
    if samplerate as sint_t > 192000 as libc::c_int * 8 as libc::c_int {
        aubio_log(AUBIO_LOG_ERR as libc::c_int,
                  b"AUBIO ERROR: %s: failed creating %s, samplerate %dHz is too large\n\x00"
                      as *const u8 as *const libc::c_char, kind, path,
                  samplerate);
        return AUBIO_FAIL as libc::c_int as uint_t
    }
    return AUBIO_OK as libc::c_int as uint_t;
}
#[no_mangle]
pub unsafe extern "C" fn aubio_io_validate_channels(mut kind: *const char_t,
                                                    mut path: *const char_t,
                                                    mut channels: uint_t)
 -> uint_t {
    if channels as sint_t <= 0 as libc::c_int {
        aubio_log(AUBIO_LOG_ERR as libc::c_int,
                  b"AUBIO ERROR: sink_%s: failed creating %s, channels should be positive, not %d\n\x00"
                      as *const u8 as *const libc::c_char, kind, path,
                  channels);
        return AUBIO_FAIL as libc::c_int as uint_t
    }
    if channels > 1024 as libc::c_int as libc::c_uint {
        aubio_log(AUBIO_LOG_ERR as libc::c_int,
                  b"AUBIO ERROR: sink_%s: failed creating %s, too many channels (%d but %d available)\n\x00"
                      as *const u8 as *const libc::c_char, kind, path,
                  channels, 1024 as libc::c_int);
        return AUBIO_FAIL as libc::c_int as uint_t
    }
    return AUBIO_OK as libc::c_int as uint_t;
}
#[no_mangle]
pub unsafe extern "C" fn aubio_source_validate_input_length(mut kind:
                                                                *const char_t,
                                                            mut path:
                                                                *const char_t,
                                                            mut hop_size:
                                                                uint_t,
                                                            mut read_data_length:
                                                                uint_t)
 -> uint_t {
    let mut length: uint_t = hop_size;
    if hop_size < read_data_length {
        aubio_log(AUBIO_LOG_WRN as libc::c_int,
                  b"AUBIO WARNING: %s: partial read from %s, trying to read %d frames, but hop_size is %d\n\x00"
                      as *const u8 as *const libc::c_char, kind, path,
                  read_data_length, hop_size);
    } else if hop_size > read_data_length {
        aubio_log(AUBIO_LOG_WRN as libc::c_int,
                  b"AUBIO WARNING: %s: partial read from %s, trying to read %d frames into a buffer of length %d\n\x00"
                      as *const u8 as *const libc::c_char, kind, path,
                  hop_size, read_data_length);
        length = read_data_length
    }
    return length;
}
#[no_mangle]
pub unsafe extern "C" fn aubio_source_validate_input_channels(mut kind:
                                                                  *const char_t,
                                                              mut path:
                                                                  *const char_t,
                                                              mut source_channels:
                                                                  uint_t,
                                                              mut read_data_height:
                                                                  uint_t)
 -> uint_t {
    let mut channels: uint_t = source_channels;
    if read_data_height < source_channels {
        aubio_log(AUBIO_LOG_WRN as libc::c_int,
                  b"AUBIO WARNING: %s: partial read from %s, trying to read %d channels, but found output of height %d\n\x00"
                      as *const u8 as *const libc::c_char, kind, path,
                  source_channels, read_data_height);
        channels = read_data_height
    } else if read_data_height > source_channels {
        // do not show a warning when trying to read into more channels than
    // the input source.
        channels = source_channels
    }
    return channels;
}
#[no_mangle]
pub unsafe extern "C" fn aubio_source_pad_output(mut read_data: *mut fvec_t,
                                                 mut source_read: uint_t) {
    if source_read < (*read_data).length {
        memset((*read_data).data.offset(source_read as isize) as
                   *mut libc::c_void, 0 as libc::c_int,
               ((*read_data).length.wrapping_sub(source_read) as
                    libc::c_ulong).wrapping_mul(::std::mem::size_of::<smpl_t>()
                                                    as libc::c_ulong));
    };
}
#[no_mangle]
pub unsafe extern "C" fn aubio_source_pad_multi_output(mut read_data:
                                                           *mut fmat_t,
                                                       mut source_channels:
                                                           uint_t,
                                                       mut source_read:
                                                           uint_t) {
    let mut i: uint_t = 0;
    if source_read < (*read_data).length {
        i = 0 as libc::c_int as uint_t;
        while i < (*read_data).height {
            memset((*(*read_data).data.offset(i as
                                                  isize)).offset(source_read
                                                                     as isize)
                       as *mut libc::c_void, 0 as libc::c_int,
                   ((*read_data).length.wrapping_sub(source_read) as
                        libc::c_ulong).wrapping_mul(::std::mem::size_of::<smpl_t>()
                                                        as libc::c_ulong));
            i = i.wrapping_add(1)
        }
    }
    // destination matrix has more channels than the file
  // copy channels from the source to extra output channels
    if (*read_data).height > source_channels {
        i = source_channels;
        while i < (*read_data).height {
            memcpy(*(*read_data).data.offset(i as isize) as *mut libc::c_void,
                   *(*read_data).data.offset(i.wrapping_rem(source_channels)
                                                 as isize) as
                       *const libc::c_void,
                   (::std::mem::size_of::<smpl_t>() as
                        libc::c_ulong).wrapping_mul((*read_data).length as
                                                        libc::c_ulong));
            i = i.wrapping_add(1)
        }
    };
}
#[no_mangle]
pub unsafe extern "C" fn aubio_sink_validate_input_length(mut kind:
                                                              *const char_t,
                                                          mut path:
                                                              *const char_t,
                                                          mut max_size:
                                                              uint_t,
                                                          mut write_data_length:
                                                              uint_t,
                                                          mut write: uint_t)
 -> uint_t {
    let mut can_write: uint_t = write;
    if write > max_size {
        aubio_log(AUBIO_LOG_WRN as libc::c_int,
                  b"AUBIO WARNING: %s: partial write to %s, trying to write %d frames, at most %d can be written at once\n\x00"
                      as *const u8 as *const libc::c_char, kind, path, write,
                  max_size);
        can_write = max_size
    }
    if can_write > write_data_length {
        aubio_log(AUBIO_LOG_WRN as libc::c_int,
                  b"AUBIO WARNING: %s: partial write to %s, trying to write %d frames, but found input of length %d\n\x00"
                      as *const u8 as *const libc::c_char, kind, path, write,
                  write_data_length);
        can_write = write_data_length
    }
    return can_write;
}
#[no_mangle]
pub unsafe extern "C" fn aubio_sink_validate_input_channels(mut kind:
                                                                *const char_t,
                                                            mut path:
                                                                *const char_t,
                                                            mut sink_channels:
                                                                uint_t,
                                                            mut write_data_height:
                                                                uint_t)
 -> uint_t {
    let mut channels: uint_t = sink_channels;
    if write_data_height < sink_channels {
        aubio_log(AUBIO_LOG_WRN as libc::c_int,
                  b"AUBIO WARNING: %s: partial write to %s, trying to write %d channels, but found input of height %d\n\x00"
                      as *const u8 as *const libc::c_char, kind, path,
                  sink_channels, write_data_height);
        channels = write_data_height
    }
    return channels;
}
