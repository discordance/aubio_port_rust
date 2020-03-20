use ::libc;
extern "C" {
    pub type _aubio_pitch_t;
    pub type _aubio_onset_t;
    #[no_mangle]
    fn calloc(_: libc::c_ulong, _: libc::c_ulong) -> *mut libc::c_void;
    #[no_mangle]
    fn free(_: *mut libc::c_void);
    #[no_mangle]
    fn aubio_log(level: sint_t, fmt: *const char_t, _: ...) -> uint_t;
    #[no_mangle]
    fn fvec_median(v: *mut fvec_t) -> smpl_t;
    #[no_mangle]
    fn aubio_level_detection(v: *const fvec_t, threshold: smpl_t) -> smpl_t;
    #[no_mangle]
    fn fvec_copy(s: *const fvec_t, t: *mut fvec_t);
    #[no_mangle]
    fn fvec_zeros(s: *mut fvec_t);
    #[no_mangle]
    fn del_fvec(s: *mut fvec_t);
    #[no_mangle]
    fn new_fvec(length: uint_t) -> *mut fvec_t;
    #[no_mangle]
    fn floorf(_: libc::c_float) -> libc::c_float;
    #[no_mangle]
    fn floor(_: libc::c_double) -> libc::c_double;
    #[no_mangle]
    fn strcmp(_: *const libc::c_char, _: *const libc::c_char) -> libc::c_int;
    /* * execute pitch detection on an input signal frame

  \param o pitch detection object as returned by new_aubio_pitch()
  \param in input signal of size [hop_size]
  \param out output pitch candidates of size [1]

*/
    #[no_mangle]
    fn aubio_pitch_do(o: *mut aubio_pitch_t, in_0: *const fvec_t,
                      out: *mut fvec_t);
    /* * change yin or yinfft tolerance threshold

  \param o pitch detection object as returned by new_aubio_pitch()
  \param tol tolerance default is 0.15 for yin and 0.85 for yinfft

*/
    #[no_mangle]
    fn aubio_pitch_set_tolerance(o: *mut aubio_pitch_t, tol: smpl_t)
     -> uint_t;
    /* * deletion of the pitch detection object

  \param o pitch detection object as returned by new_aubio_pitch()

*/
    #[no_mangle]
    fn del_aubio_pitch(o: *mut aubio_pitch_t);
    /* * creation of the pitch detection object

  \param method set pitch detection algorithm
  \param buf_size size of the input buffer to analyse
  \param hop_size step size between two consecutive analysis instant
  \param samplerate sampling rate of the signal

  \return newly created ::aubio_pitch_t

*/
    #[no_mangle]
    fn new_aubio_pitch(method: *const char_t, buf_size: uint_t,
                       hop_size: uint_t, samplerate: uint_t)
     -> *mut aubio_pitch_t;
    /* * set the output unit of the pitch detection object

  \param o pitch detection object as returned by new_aubio_pitch()
  \param mode set pitch units for output

  mode can be one of "Hz", "midi", "cent", or "bin". Defaults to "Hz".

  \return 0 if successfull, non-zero otherwise

*/
    #[no_mangle]
    fn aubio_pitch_set_unit(o: *mut aubio_pitch_t, mode: *const char_t)
     -> uint_t;
    /* * set the silence threshold of the pitch detection object

  \param o pitch detection object as returned by new_aubio_pitch()
  \param silence level threshold under which pitch should be ignored, in dB

  \return 0 if successfull, non-zero otherwise

*/
    #[no_mangle]
    fn aubio_pitch_set_silence(o: *mut aubio_pitch_t, silence: smpl_t)
     -> uint_t;
    /* * set the silence threshold of the pitch detection object

  \param o pitch detection object as returned by ::new_aubio_pitch()

  \return level threshold under which pitch should be ignored, in dB

*/
    #[no_mangle]
    fn aubio_pitch_get_silence(o: *mut aubio_pitch_t) -> smpl_t;
    /* * create onset detection object

  \param method onset detection type as specified in specdesc.h
  \param buf_size buffer size for phase vocoder
  \param hop_size hop size for phase vocoder
  \param samplerate sampling rate of the input signal

  \return newly created ::aubio_onset_t

*/
    #[no_mangle]
    fn new_aubio_onset(method: *const char_t, buf_size: uint_t,
                       hop_size: uint_t, samplerate: uint_t)
     -> *mut aubio_onset_t;
    /* * execute onset detection

  \param o onset detection object as returned by new_aubio_onset()
  \param input new audio vector of length hop_size
  \param onset output vector of length 1, containing 0 if no onset was found,
  and a value equal or greater than 1 otherwise

  When no onset was detected, the first element of the output vector `onset`
  is set to 0.

  When an onset is found, the first element of the output vector `onset` is set
  to `offset = 1 + a` where `a` is a number in the range`[0, 1]`.

  The final onset detection time, in samples, can be obtained with
  aubio_onset_get_last(). It can also be derived from `offset` as
  follows:

  \code
    t = total_frames + offset * hop_size - delay
  \endcode

  where `total_frames` is the total number of frames processed so far, and
  `delay` is the current delay of the onset object, as returned by
  aubio_onset_get_delay().

*/
    #[no_mangle]
    fn aubio_onset_do(o: *mut aubio_onset_t, input: *const fvec_t,
                      onset: *mut fvec_t);
    /* * set onset detection silence threshold

  \param o onset detection object as returned by new_aubio_onset()
  \param silence new silence detection threshold

*/
    #[no_mangle]
    fn aubio_onset_set_silence(o: *mut aubio_onset_t, silence: smpl_t)
     -> uint_t;
    /* * set onset detection peak picking threshold

  \param o onset detection object as returned by new_aubio_onset()
  \param threshold new peak-picking threshold

*/
    #[no_mangle]
    fn aubio_onset_set_threshold(o: *mut aubio_onset_t, threshold: smpl_t)
     -> uint_t;
    /* * set minimum inter onset interval in milliseconds

  \param o onset detection object as returned by new_aubio_onset()
  \param minioi minimum interval between two consecutive onsets (in
  milliseconds)

*/
    #[no_mangle]
    fn aubio_onset_set_minioi_ms(o: *mut aubio_onset_t, minioi: smpl_t)
     -> uint_t;
    /* * get minimum inter onset interval in milliseconds

  \param o onset detection object as returned by new_aubio_onset()
  \return minimum interval between two consecutive onsets (in
  milliseconds)

*/
    #[no_mangle]
    fn aubio_onset_get_minioi_ms(o: *const aubio_onset_t) -> smpl_t;
    /* * delete onset detection object

  \param o onset detection object to delete

*/
    #[no_mangle]
    fn del_aubio_onset(o: *mut aubio_onset_t);
}
pub type smpl_t = libc::c_float;
pub type uint_t = libc::c_uint;
pub type sint_t = libc::c_int;
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
/* * print format for sample in single precision */
/* * long sample format (64 bits or more) */
/* * print format for sample in double precision */
/* * unsigned integer */
/* * signed integer */
/* * character */
pub type char_t = libc::c_char;
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
/*
  Copyright (C) 2003-2013 Paul Brossier <piem@aubio.org>

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

  Pitch detection object

  This file creates the objects required for the computation of the selected
  pitch detection algorithm and output the results, in midi note or Hz.

  \section pitch Pitch detection methods

  A list of the pitch detection methods currently available follows.

  \b \p default : use the default method

  Currently, the default method is set to \p yinfft .

  \b \p schmitt : Schmitt trigger

  This pitch extraction method implements a Schmitt trigger to estimate the
  period of a signal.

  This file was derived from the tuneit project, written by Mario Lang to
  detect the fundamental frequency of a sound.

  See http://delysid.org/tuneit.html

  \b \p fcomb : a fast harmonic comb filter

  This pitch extraction method implements a fast harmonic comb filter to
  determine the fundamental frequency of a harmonic sound.

  This file was derived from the tuneit project, written by Mario Lang to
  detect the fundamental frequency of a sound.

  See http://delysid.org/tuneit.html

  \b \p mcomb : multiple-comb filter

  This fundamental frequency estimation algorithm implements spectral
  flattening, multi-comb filtering and peak histogramming.

  This method was designed by Juan P. Bello and described in:

  Juan-Pablo Bello. ``Towards the Automated Analysis of Simple Polyphonic
  Music''.  PhD thesis, Centre for Digital Music, Queen Mary University of
  London, London, UK, 2003.

  \b \p yin : YIN algorithm

  This algorithm was developed by A. de Cheveigne and H. Kawahara and
  published in:

  De Cheveigné, A., Kawahara, H. (2002) "YIN, a fundamental frequency
  estimator for speech and music", J. Acoust. Soc. Am. 111, 1917-1930.

  see http://recherche.ircam.fr/equipes/pcm/pub/people/cheveign.html

  \b \p yinfast : Yinfast algorithm

  This algorithm is equivalent to the YIN algorithm, but computed in the
  spectral domain for efficiency. See also `python/demos/demo_yin_compare.py`.

  \b \p yinfft : Yinfft algorithm

  This algorithm was derived from the YIN algorithm. In this implementation, a
  Fourier transform is used to compute a tapered square difference function,
  which allows spectral weighting. Because the difference function is tapered,
  the selection of the period is simplified.

  Paul Brossier, [Automatic annotation of musical audio for interactive
  systems](http://aubio.org/phd/), Chapter 3, Pitch Analysis, PhD thesis,
  Centre for Digital music, Queen Mary University of London, London, UK, 2006.

  \example pitch/test-pitch.c
  \example examples/aubiopitch.c

*/
/* * pitch detection object */
pub type aubio_pitch_t = _aubio_pitch_t;
/*
  Copyright (C) 2006-2013 Paul Brossier <piem@aubio.org>

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

  Onset detection object

  The following routines compute the onset detection function and detect peaks
  in these functions. When onsets are found above a given silence threshold,
  and after a minimum inter-onset interval, the output vector returned by
  aubio_onset_do() is filled with `1`. Otherwise, the output vector remains
  `0`.

  The peak-picking threshold, the silence threshold, and the minimum
  inter-onset interval can be adjusted during the execution of the
  aubio_onset_do routine using the corresponding functions.

  \example onset/test-onset.c
  \example examples/aubioonset.c
  \example examples/aubionotes.c

*/
/* * onset detection object */
pub type aubio_onset_t = _aubio_onset_t;
#[derive(Copy, Clone)]
#[repr(C)]
pub struct _aubio_notes_t {
    pub onset_buf_size: uint_t,
    pub pitch_buf_size: uint_t,
    pub hop_size: uint_t,
    pub samplerate: uint_t,
    pub median: uint_t,
    pub note_buffer: *mut fvec_t,
    pub note_buffer2: *mut fvec_t,
    pub pitch: *mut aubio_pitch_t,
    pub pitch_output: *mut fvec_t,
    pub pitch_tolerance: smpl_t,
    pub onset: *mut aubio_onset_t,
    pub onset_output: *mut fvec_t,
    pub onset_threshold: smpl_t,
    pub curnote: smpl_t,
    pub newnote: smpl_t,
    pub silence_threshold: smpl_t,
    pub isready: uint_t,
    pub last_onset_level: smpl_t,
    pub release_drop_level: smpl_t,
}
/*
  Copyright (C) 2003-2014 Paul Brossier <piem@aubio.org>

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

  Note detection object

*/
/* * notes detection object */
pub type aubio_notes_t = _aubio_notes_t;
/* * create notes detection object

  \param method notes detection type as specified in specdesc.h
  \param buf_size buffer size for phase vocoder
  \param hop_size hop size for phase vocoder
  \param samplerate sampling rate of the input signal

  \return newly created ::aubio_notes_t

*/
#[no_mangle]
pub unsafe extern "C" fn new_aubio_notes(mut method: *const char_t,
                                         mut buf_size: uint_t,
                                         mut hop_size: uint_t,
                                         mut samplerate: uint_t)
 -> *mut aubio_notes_t {
    let mut o: *mut aubio_notes_t =
        calloc(::std::mem::size_of::<aubio_notes_t>() as libc::c_ulong,
               1 as libc::c_int as libc::c_ulong) as *mut aubio_notes_t;
    let mut onset_method: *const char_t =
        b"default\x00" as *const u8 as *const libc::c_char;
    let mut pitch_method: *const char_t =
        b"default\x00" as *const u8 as *const libc::c_char;
    (*o).onset_buf_size = buf_size;
    (*o).pitch_buf_size =
        buf_size.wrapping_mul(4 as libc::c_int as libc::c_uint);
    (*o).hop_size = hop_size;
    (*o).onset_threshold = 0.0f64 as smpl_t;
    (*o).pitch_tolerance = 0.0f64 as smpl_t;
    (*o).samplerate = samplerate;
    (*o).median = 6 as libc::c_int as uint_t;
    (*o).isready = 0 as libc::c_int as uint_t;
    (*o).onset =
        new_aubio_onset(onset_method, (*o).onset_buf_size, (*o).hop_size,
                        (*o).samplerate);
    if !(*o).onset.is_null() {
        if (*o).onset_threshold as libc::c_double != 0.0f64 {
            aubio_onset_set_threshold((*o).onset, (*o).onset_threshold);
        }
        (*o).onset_output = new_fvec(1 as libc::c_int as uint_t);
        (*o).pitch =
            new_aubio_pitch(pitch_method, (*o).pitch_buf_size, (*o).hop_size,
                            (*o).samplerate);
        if !(*o).pitch.is_null() {
            if (*o).pitch_tolerance as libc::c_double != 0.0f64 {
                aubio_pitch_set_tolerance((*o).pitch, (*o).pitch_tolerance);
            }
            aubio_pitch_set_unit((*o).pitch,
                                 b"midi\x00" as *const u8 as
                                     *const libc::c_char);
            (*o).pitch_output = new_fvec(1 as libc::c_int as uint_t);
            if strcmp(method,
                      b"default\x00" as *const u8 as *const libc::c_char) !=
                   0 as libc::c_int {
                aubio_log(AUBIO_LOG_ERR as libc::c_int,
                          b"AUBIO ERROR: notes: unknown notes detection method \"%s\"\n\x00"
                              as *const u8 as *const libc::c_char, method);
            } else {
                (*o).note_buffer = new_fvec((*o).median);
                (*o).note_buffer2 = new_fvec((*o).median);
                if !((*o).onset_output.is_null() ||
                         (*o).pitch_output.is_null() ||
                         (*o).note_buffer.is_null() ||
                         (*o).note_buffer2.is_null()) {
                    (*o).curnote = -1.0f64 as smpl_t;
                    (*o).newnote = 0.0f64 as smpl_t;
                    aubio_notes_set_silence(o, -70.0f64 as smpl_t);
                    aubio_notes_set_minioi_ms(o, 30.0f64 as smpl_t);
                    (*o).last_onset_level = -70.0f64 as smpl_t;
                    (*o).release_drop_level = 10.0f64 as smpl_t;
                    return o
                }
            }
        }
    }
    del_aubio_notes(o);
    return 0 as *mut aubio_notes_t;
}
/* * set notes detection silence threshold

  \param o notes detection object as returned by new_aubio_notes()
  \param silence new silence detection threshold

  \return 0 on success, non-zero otherwise

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_notes_set_silence(mut o: *mut aubio_notes_t,
                                                 mut silence: smpl_t)
 -> uint_t {
    let mut err: uint_t = AUBIO_OK as libc::c_int as uint_t;
    if aubio_pitch_set_silence((*o).pitch, silence) !=
           AUBIO_OK as libc::c_int as libc::c_uint {
        err = AUBIO_FAIL as libc::c_int as uint_t
    }
    if aubio_onset_set_silence((*o).onset, silence) !=
           AUBIO_OK as libc::c_int as libc::c_uint {
        err = AUBIO_FAIL as libc::c_int as uint_t
    }
    (*o).silence_threshold = silence;
    return err;
}
/* * get notes detection silence threshold

  \param o notes detection object as returned by new_aubio_notes()

  \return current silence threshold

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_notes_get_silence(mut o: *const aubio_notes_t)
 -> smpl_t {
    return aubio_pitch_get_silence((*o).pitch);
}
/* * set notes detection minimum inter-onset interval, in millisecond

  \param o notes detection object as returned by new_aubio_notes()
  \param minioi_ms new inter-onset interval

  \return 0 on success, non-zero otherwise

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_notes_set_minioi_ms(mut o: *mut aubio_notes_t,
                                                   mut minioi_ms: smpl_t)
 -> uint_t {
    let mut err: uint_t = AUBIO_OK as libc::c_int as uint_t;
    if (*o).onset.is_null() ||
           aubio_onset_set_minioi_ms((*o).onset, minioi_ms) !=
               0 as libc::c_int as libc::c_uint {
        err = AUBIO_FAIL as libc::c_int as uint_t
    }
    return err;
}
/* * get notes detection minimum inter-onset interval, in millisecond

  \param o notes detection object as returned by new_aubio_notes()

  \return current minimum inter onset interval

 */
#[no_mangle]
pub unsafe extern "C" fn aubio_notes_get_minioi_ms(mut o:
                                                       *const aubio_notes_t)
 -> smpl_t {
    return aubio_onset_get_minioi_ms((*o).onset);
}
/* * set note release drop level, in dB

  This function sets the release_drop_level parameter, in dB. When a new note
  is found, the current level in dB is measured. If the measured level drops
  under that initial level - release_drop_level, then a note-off will be
  emitted.

  Defaults to `10`, in dB.

  \note This parameter was added in version `0.4.8`. Results obtained with
  earlier versions can be reproduced by setting this value to `100`, so that
  note-off will not be played until the next note.

  \param o notes detection object as returned by new_aubio_notes()
  \param release_drop new release drop level, in dB

  \return 0 on success, non-zero otherwise

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_notes_set_release_drop(mut o:
                                                          *mut aubio_notes_t,
                                                      mut release_drop_level:
                                                          smpl_t) -> uint_t {
    let mut err: uint_t = AUBIO_OK as libc::c_int as uint_t;
    if release_drop_level as libc::c_double <= 0.0f64 {
        aubio_log(AUBIO_LOG_ERR as libc::c_int,
                  b"AUBIO ERROR: notes: release_drop should be >= 0, got %f\n\x00"
                      as *const u8 as *const libc::c_char,
                  release_drop_level as libc::c_double);
        err = AUBIO_FAIL as libc::c_int as uint_t
    } else { (*o).release_drop_level = release_drop_level }
    return err;
}
/* * get notes object release drop level, in dB

  \param o notes detection object as returned by new_aubio_notes()

  \return current release drop level, in dB

 */
#[no_mangle]
pub unsafe extern "C" fn aubio_notes_get_release_drop(mut o:
                                                          *const aubio_notes_t)
 -> smpl_t {
    return (*o).release_drop_level;
}
/* * append new note candidate to the note_buffer and return filtered value. we
 * need to copy the input array as fvec_median destroy its input data.*/
unsafe extern "C" fn note_append(mut note_buffer: *mut fvec_t,
                                 mut curnote: smpl_t) {
    let mut i: uint_t = 0 as libc::c_int as uint_t;
    i = 0 as libc::c_int as uint_t;
    while i <
              (*note_buffer).length.wrapping_sub(1 as libc::c_int as
                                                     libc::c_uint) {
        *(*note_buffer).data.offset(i as isize) =
            *(*note_buffer).data.offset(i.wrapping_add(1 as libc::c_int as
                                                           libc::c_uint) as
                                            isize);
        i = i.wrapping_add(1)
    }
    //note_buffer->data[note_buffer->length - 1] = ROUND(10.*curnote)/10.;
    *(*note_buffer).data.offset((*note_buffer).length.wrapping_sub(1 as
                                                                       libc::c_int
                                                                       as
                                                                       libc::c_uint)
                                    as isize) =
        floorf((1.0f64 * curnote as libc::c_double + 0.5f64) as
                   libc::c_float);
}
unsafe extern "C" fn aubio_notes_get_latest_note(mut o: *mut aubio_notes_t)
 -> smpl_t {
    fvec_copy((*o).note_buffer, (*o).note_buffer2);
    return (fvec_median((*o).note_buffer2) as libc::c_double / 1.0f64) as
               smpl_t;
}
/* * execute note detection on an input signal frame

  \param o note detection object as returned by new_aubio_notes()
  \param input input signal of size [hop_size]
  \param output output notes, fvec of length 3

  The notes output is a vector of length 3 containing:
   - 0. the midi note value, or 0 if no note was found
   - 1. the note velocity
   - 2. the midi note to turn off

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_notes_do(mut o: *mut aubio_notes_t,
                                        mut input: *const fvec_t,
                                        mut notes: *mut fvec_t) {
    let mut new_pitch: smpl_t = 0.;
    let mut curlevel: smpl_t = 0.;
    fvec_zeros(notes);
    aubio_onset_do((*o).onset, input, (*o).onset_output);
    aubio_pitch_do((*o).pitch, input, (*o).pitch_output);
    new_pitch = *(*(*o).pitch_output).data.offset(0 as libc::c_int as isize);
    if (*o).median != 0 { note_append((*o).note_buffer, new_pitch); }
    /* curlevel is negatif or 1 if silence */
    curlevel = aubio_level_detection(input, (*o).silence_threshold);
    if *(*(*o).onset_output).data.offset(0 as libc::c_int as isize) !=
           0 as libc::c_int as libc::c_float {
        /* test for silence */
        if curlevel as libc::c_double == 1.0f64 {
            if (*o).median != 0 { (*o).isready = 0 as libc::c_int as uint_t }
            /* send note off */
      //send_noteon(o->curnote,0);
      //notes->data[0] = o->curnote;
      //notes->data[1] = 0.;
      //AUBIO_WRN("notes: sending note-off at onset, not enough level\n");
            *(*notes).data.offset(2 as libc::c_int as isize) = (*o).curnote
        } else {
            if (*o).median != 0 {
                (*o).isready = 1 as libc::c_int as uint_t
            } else {
                /* kill old note */
        //send_noteon(o->curnote,0, o->samplerate);
        //AUBIO_WRN("notes: sending note-off at onset, new onset detected\n");
                *(*notes).data.offset(2 as libc::c_int as isize) =
                    (*o).curnote;
                /* get and send new one */
        //send_noteon(new_pitch,127+(int)floor(curlevel), o->samplerate);
                *(*notes).data.offset(0 as libc::c_int as isize) = new_pitch;
                *(*notes).data.offset(1 as libc::c_int as isize) =
                    (127 as libc::c_int +
                         floor(curlevel as libc::c_double) as libc::c_int) as
                        smpl_t;
                (*o).curnote = new_pitch
            }
            (*o).last_onset_level = curlevel
        }
    } else if curlevel < (*o).last_onset_level - (*o).release_drop_level {
        // send note off
      //AUBIO_WRN("notes: sending note-off, release detected\n");
        *(*notes).data.offset(0 as libc::c_int as isize) =
            0 as libc::c_int as smpl_t;
        *(*notes).data.offset(1 as libc::c_int as isize) =
            0 as libc::c_int as smpl_t;
        *(*notes).data.offset(2 as libc::c_int as isize) = (*o).curnote;
        // reset last_onset_level to silence_threshold
        (*o).last_onset_level = (*o).silence_threshold;
        (*o).curnote = 0 as libc::c_int as smpl_t
    } else if (*o).median != 0 {
        if (*o).isready > 0 as libc::c_int as libc::c_uint {
            (*o).isready = (*o).isready.wrapping_add(1)
        }
        if (*o).isready == (*o).median {
            /* kill old note */
        //send_noteon(curnote,0);
            if (*o).curnote != 0 as libc::c_int as libc::c_float {
                //AUBIO_WRN("notes: sending note-off, new note detected\n");
                *(*notes).data.offset(2 as libc::c_int as isize) =
                    (*o).curnote
            }
            (*o).newnote = aubio_notes_get_latest_note(o);
            (*o).curnote = (*o).newnote;
            /* get and send new one */
            if (*o).curnote > 45 as libc::c_int as libc::c_float {
                //send_noteon(curnote,127+(int)floor(curlevel));
                *(*notes).data.offset(0 as libc::c_int as isize) =
                    (*o).curnote;
                *(*notes).data.offset(1 as libc::c_int as isize) =
                    (127 as libc::c_int +
                         floor(curlevel as libc::c_double) as libc::c_int) as
                        smpl_t
            }
        }
    };
}
// if median
/* * delete notes detection object

  \param o notes detection object to delete

*/
#[no_mangle]
pub unsafe extern "C" fn del_aubio_notes(mut o: *mut aubio_notes_t) {
    if !(*o).note_buffer.is_null() { del_fvec((*o).note_buffer); }
    if !(*o).note_buffer2.is_null() { del_fvec((*o).note_buffer2); }
    if !(*o).pitch_output.is_null() { del_fvec((*o).pitch_output); }
    if !(*o).pitch.is_null() { del_aubio_pitch((*o).pitch); }
    if !(*o).onset_output.is_null() { del_fvec((*o).onset_output); }
    if !(*o).onset.is_null() { del_aubio_onset((*o).onset); }
    free(o as *mut libc::c_void);
}
