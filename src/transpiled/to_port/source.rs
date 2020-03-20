use ::libc;
extern "C" {
    #[no_mangle]
    fn calloc(_: libc::c_ulong, _: libc::c_ulong) -> *mut libc::c_void;
    #[no_mangle]
    fn free(_: *mut libc::c_void);
    /* file interface */
    /* strings */
    /* Error reporting */
    /* Logging */
    /* * internal logging function, defined in utils/log.c */
    #[no_mangle]
    fn aubio_log(level: sint_t, fmt: *const char_t, _: ...) -> uint_t;
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
#[derive(Copy, Clone)]
#[repr(C)]
pub struct _aubio_source_t {
    pub source: *mut libc::c_void,
    pub s_do: aubio_source_do_t,
    pub s_do_multi: aubio_source_do_multi_t,
    pub s_get_samplerate: aubio_source_get_samplerate_t,
    pub s_get_channels: aubio_source_get_channels_t,
    pub s_get_duration: aubio_source_get_duration_t,
    pub s_seek: aubio_source_seek_t,
    pub s_close: aubio_source_close_t,
    pub s_del: del_aubio_source_t,
}
pub type del_aubio_source_t
    =
    Option<unsafe extern "C" fn(_: *mut aubio_source_t) -> ()>;
/*
  Copyright (C) 2012-2013 Paul Brossier <piem@aubio.org>

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

  Media source to read blocks of consecutive audio samples from file.

  To write to file, use ::aubio_sink_t.

  Depending on how aubio was compiled, the following sources will be available.

  When creating a new source using ::new_aubio_source, the new function of each
  of the compiled-in sources will be used, in the following order, until one of
  them gets successfully created. If all sources returned NULL,
  ::new_aubio_source will return NULL.

  \b \p source_avcodec : libav

  aubio can be optionally compiled with [libav](http://libav.org), which can
  read from a very large number of audio and video formats, including over
  different network protocols such as HTTP.

  \b \p source_apple_audio : ExtAudioFileRef

  On Mac and iOS platforms, aubio should be compiled with CoreAudio [Extended
  Audio File Services]
  (https://developer.apple.com/library/mac/documentation/MusicAudio/Reference/ExtendedAudioFileServicesReference/Reference/reference.html).
  This provides access to most common audio file formats, including compressed
  ones.

  \b \p source_sndfile : libsndfile

  Also optional, aubio can be built against
  [libsndfile](http://www.mega-nerd.com/libsndfile/), which can read [most
  uncompressed formats](http://www.mega-nerd.com/libsndfile/#Features).

  \b \p source_wavread : native WAV reader

  A simple source to read from 16-bits PCM RIFF encoded WAV files.

  \example io/test-source.c

*/
/* * media source object */
pub type aubio_source_t = _aubio_source_t;
pub type aubio_source_close_t
    =
    Option<unsafe extern "C" fn(_: *mut aubio_source_t) -> uint_t>;
pub type aubio_source_seek_t
    =
    Option<unsafe extern "C" fn(_: *mut aubio_source_t, _: uint_t) -> uint_t>;
pub type aubio_source_get_duration_t
    =
    Option<unsafe extern "C" fn(_: *mut aubio_source_t) -> uint_t>;
pub type aubio_source_get_channels_t
    =
    Option<unsafe extern "C" fn(_: *mut aubio_source_t) -> uint_t>;
pub type aubio_source_get_samplerate_t
    =
    Option<unsafe extern "C" fn(_: *mut aubio_source_t) -> uint_t>;
pub type aubio_source_do_multi_t
    =
    Option<unsafe extern "C" fn(_: *mut aubio_source_t, _: *mut fmat_t,
                                _: *mut uint_t) -> ()>;
/*
  Copyright (C) 2012 Paul Brossier <piem@aubio.org>

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
/* HAVE_LIBAV */
/* HAVE_SOURCE_APPLE_AUDIO */
/* HAVE_SNDFILE */
/* HAVE_WAVREAD */
pub type aubio_source_do_t
    =
    Option<unsafe extern "C" fn(_: *mut aubio_source_t, _: *mut fvec_t,
                                _: *mut uint_t) -> ()>;
/* *

  create new ::aubio_source_t

  \param uri the file path or uri to read from
  \param samplerate sampling rate to view the fie at
  \param hop_size the size of the blocks to read from

  Creates a new source object. If `0` is passed as `samplerate`, the sample
  rate of the original file is used.

  The samplerate of newly created source can be obtained using
  ::aubio_source_get_samplerate.

*/
#[no_mangle]
pub unsafe extern "C" fn new_aubio_source(mut uri: *const char_t,
                                          mut samplerate: uint_t,
                                          mut hop_size: uint_t)
 -> *mut aubio_source_t {
    let mut s: *mut aubio_source_t =
        calloc(::std::mem::size_of::<aubio_source_t>() as libc::c_ulong,
               1 as libc::c_int as libc::c_ulong) as *mut aubio_source_t;
    /* HAVE_LIBAV */
    /* HAVE_SOURCE_APPLE_AUDIO */
    /* HAVE_SNDFILE */
    /* HAVE_WAVREAD */
    aubio_log(AUBIO_LOG_ERR as libc::c_int,
              b"AUBIO ERROR: source: failed creating with %s at %dHz with hop size %d (no source built-in)\n\x00"
                  as *const u8 as *const libc::c_char, uri, samplerate,
              hop_size);
    del_aubio_source(s);
    return 0 as *mut aubio_source_t;
}
/* *

  read monophonic vector of length hop_size from source object

  \param s source object, created with ::new_aubio_source
  \param read_to ::fvec_t of data to read to
  \param read upon returns, equals to number of frames actually read

  Upon returns, `read` contains the number of frames actually read from the
  source. `hop_size` if enough frames could be read, less otherwise.

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_source_do(mut s: *mut aubio_source_t,
                                         mut data: *mut fvec_t,
                                         mut read: *mut uint_t) {
    (*s).s_do.expect("non-null function pointer")((*s).source as
                                                      *mut aubio_source_t,
                                                  data, read);
}
/* *

  read polyphonic vector of length hop_size from source object

  \param s source object, created with ::new_aubio_source
  \param read_to ::fmat_t of data to read to
  \param[out] read upon returns, equals to number of frames actually read

  Upon returns, `read` contains the number of frames actually read from the
  source. `hop_size` if enough frames could be read, less otherwise.

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_source_do_multi(mut s: *mut aubio_source_t,
                                               mut data: *mut fmat_t,
                                               mut read: *mut uint_t) {
    (*s).s_do_multi.expect("non-null function pointer")((*s).source as
                                                            *mut aubio_source_t,
                                                        data, read);
}
/* *

  close source object

  \param s source object, created with ::new_aubio_source

  \return 0 if sucessful, non-zero on failure

 */
#[no_mangle]
pub unsafe extern "C" fn aubio_source_close(mut s: *mut aubio_source_t)
 -> uint_t {
    return (*s).s_close.expect("non-null function pointer")((*s).source as
                                                                *mut aubio_source_t);
}
/* *

  close source and cleanup memory

  \param s source object, created with ::new_aubio_source

*/
#[no_mangle]
pub unsafe extern "C" fn del_aubio_source(mut s: *mut aubio_source_t) {
    //AUBIO_ASSERT(s);
    if !s.is_null() && (*s).s_del.is_some() && !(*s).source.is_null() {
        (*s).s_del.expect("non-null function pointer")((*s).source as
                                                           *mut aubio_source_t);
    }
    free(s as *mut libc::c_void);
}
/* *

  get samplerate of source object

  \param s source object, created with ::new_aubio_source
  \return samplerate, in Hz

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_source_get_samplerate(mut s:
                                                         *mut aubio_source_t)
 -> uint_t {
    return (*s).s_get_samplerate.expect("non-null function pointer")((*s).source
                                                                         as
                                                                         *mut aubio_source_t);
}
/* *

  get channels of source object

  \param s source object, created with ::new_aubio_source
  \return channels

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_source_get_channels(mut s: *mut aubio_source_t)
 -> uint_t {
    return (*s).s_get_channels.expect("non-null function pointer")((*s).source
                                                                       as
                                                                       *mut aubio_source_t);
}
/* *

  get the duration of source object, in frames

  \param s source object, created with ::new_aubio_source
  \return number of frames in file

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_source_get_duration(mut s: *mut aubio_source_t)
 -> uint_t {
    return (*s).s_get_duration.expect("non-null function pointer")((*s).source
                                                                       as
                                                                       *mut aubio_source_t);
}
/* *

  seek source object

  \param s source object, created with ::new_aubio_source
  \param pos position to seek to, in frames

  \return 0 if sucessful, non-zero on failure

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_source_seek(mut s: *mut aubio_source_t,
                                           mut seek: uint_t) -> uint_t {
    return (*s).s_seek.expect("non-null function pointer")((*s).source as
                                                               *mut aubio_source_t,
                                                           seek);
}
