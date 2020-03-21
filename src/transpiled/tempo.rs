use crate::transpiled::peakpicker::aubio_peakpicker_t;
use crate::transpiled::phasevoc::aubio_pvoc_t;
use crate::transpiled::specdesc::aubio_specdesc_t;
use crate::transpiled::beattracking::aubio_beattracking_t;

extern "C" {
    #[no_mangle]
    fn calloc(_: u64, _: u64) -> *mut core::ffi::c_void;
    #[no_mangle]
    fn free(_: *mut core::ffi::c_void);
    #[no_mangle]
    fn floorf(_: f32) -> f32;
    #[no_mangle]
    fn strcmp(_: *const i8, _: *const i8) -> i32;
    #[no_mangle]
    fn strncpy(_: *mut i8, _: *const i8, _: u64)
     -> *mut i8;
    /* * fvec_t buffer creation function

  \param length the length of the buffer to create

*/
    #[no_mangle]
    fn new_fvec(length: uint_t) -> *mut fvec_t;
    /* * fvec_t buffer deletion function

  \param s buffer to delete as returned by new_fvec()

*/
    #[no_mangle]
    fn del_fvec(s: *mut fvec_t);
    /* * check if buffer level in dB SPL is under a given threshold

  \param v vector to get level from
  \param threshold threshold in dB SPL

  \return 0 if level is under the given threshold, 1 otherwise

*/
    #[no_mangle]
    fn aubio_silence_detection(v: *const fvec_t, threshold: smpl_t) -> uint_t;
    /* * return the next power of power of 2 greater than a */
    #[no_mangle]
    fn aubio_next_power_of_two(a: uint_t) -> uint_t;
    #[no_mangle]
    fn aubio_log(level: sint_t, fmt: *const char_t, _: ...) -> uint_t;
    /* * cvec_t buffer creation function

  This function creates a cvec_t structure holding two arrays of size
  [length/2+1], corresponding to the norm and phase values of the
  spectral frame. The length stored in the structure is the actual size of both
  arrays, not the length of the complex and symmetrical vector, specified as
  creation argument.

  \param length the length of the buffer to create

*/
    #[no_mangle]
    fn new_cvec(length: uint_t) -> *mut cvec_t;
    /* * cvec_t buffer deletion function

  \param s buffer to delete as returned by new_cvec()

*/
    #[no_mangle]
    fn del_cvec(s: *mut cvec_t);
    /* * execute spectral description function on a spectral frame

  Generic function to compute spectral description.

  \param o spectral description object as returned by new_aubio_specdesc()
  \param fftgrain input signal spectrum as computed by aubio_pvoc_do
  \param desc output vector (one sample long, to send to the peak picking)

*/
    #[no_mangle]
    fn aubio_specdesc_do(o: *mut aubio_specdesc_t, fftgrain: *const cvec_t,
                         desc: *mut fvec_t);
    /* * creation of a spectral description object

  \param method spectral description method
  \param buf_size length of the input spectrum frame

  The parameter \p method is a string that can be any of:

    - onset novelty functions: `complex`, `energy`, `hfc`, `kl`, `mkl`,
    `phase`, `specdiff`, `specflux`, `wphase`,

    - spectral descriptors: `centroid`, `decrease`, `kurtosis`, `rolloff`,
    `skewness`, `slope`, `spread`.

*/
    #[no_mangle]
    fn new_aubio_specdesc(method: *const char_t, buf_size: uint_t)
     -> *mut aubio_specdesc_t;
    /* * deletion of a spectral descriptor

  \param o spectral descriptor object as returned by new_aubio_specdesc()

*/
    #[no_mangle]
    fn del_aubio_specdesc(o: *mut aubio_specdesc_t);
    /* * create beat tracking object

  \param winlen length of the onset detection window
  \param hop_size number of onset detection samples [512]
  \param samplerate samplerate of the input signal

*/
    #[no_mangle]
    fn new_aubio_beattracking(winlen: uint_t, hop_size: uint_t,
                              samplerate: uint_t)
     -> *mut aubio_beattracking_t;
    /* * track the beat

  \param bt beat tracking object
  \param dfframes current input detection function frame, smoothed by
  adaptive median threshold.
  \param out stored detected beat locations

*/
    #[no_mangle]
    fn aubio_beattracking_do(bt: *mut aubio_beattracking_t,
                             dfframes: *const fvec_t, out: *mut fvec_t);
    /* * get current beat period in samples

  \param bt beat tracking object

  Returns the currently observed period, in samples, or 0 if no consistent
  value is found.

*/
    #[no_mangle]
    fn aubio_beattracking_get_period(bt: *const aubio_beattracking_t)
     -> smpl_t;
    /* * get current beat period in seconds

  \param bt beat tracking object

  Returns the currently observed period, in seconds, or 0 if no consistent
  value is found.

*/
    #[no_mangle]
    fn aubio_beattracking_get_period_s(bt: *const aubio_beattracking_t)
     -> smpl_t;
    /* * get current tempo in bpm

  \param bt beat tracking object

  Returns the currently observed tempo, in beats per minutes, or 0 if no
  consistent value is found.

*/
    #[no_mangle]
    fn aubio_beattracking_get_bpm(bt: *const aubio_beattracking_t) -> smpl_t;
    /* * get current tempo confidence

  \param bt beat tracking object

  Returns the confidence with which the tempo has been observed, 0 if no
  consistent value is found.

*/
    #[no_mangle]
    fn aubio_beattracking_get_confidence(bt: *const aubio_beattracking_t)
     -> smpl_t;
    /* * delete beat tracking object

  \param p beat tracking object

*/
    #[no_mangle]
    fn del_aubio_beattracking(p: *mut aubio_beattracking_t);
    /* * create phase vocoder object

  \param win_s size of analysis buffer (and length the FFT transform)
  \param hop_s step size between two consecutive analysis

*/
    #[no_mangle]
    fn new_aubio_pvoc(win_s: uint_t, hop_s: uint_t) -> *mut aubio_pvoc_t;
    /* * delete phase vocoder object

  \param pv phase vocoder object as returned by new_aubio_pvoc

*/
    #[no_mangle]
    fn del_aubio_pvoc(pv: *mut aubio_pvoc_t);
    /* * compute spectral frame

  This function accepts an input vector of size [hop_s]. The
  analysis buffer is rotated and filled with the new data. After windowing of
  this signal window, the Fourier transform is computed and returned in
  fftgrain as two vectors, magnitude and phase.

  \param pv phase vocoder object as returned by new_aubio_pvoc
  \param in new input signal (hop_s long)
  \param fftgrain output spectral frame

*/
    #[no_mangle]
    fn aubio_pvoc_do(pv: *mut aubio_pvoc_t, in_0: *const fvec_t,
                     fftgrain: *mut cvec_t);
    /* * peak-picker creation function */
    #[no_mangle]
    fn new_aubio_peakpicker() -> *mut aubio_peakpicker_t;
    /* * real time peak picking function */
    #[no_mangle]
    fn aubio_peakpicker_do(p: *mut aubio_peakpicker_t, in_0: *mut fvec_t,
                           out: *mut fvec_t);
    /* * destroy peak picker structure */
    #[no_mangle]
    fn del_aubio_peakpicker(p: *mut aubio_peakpicker_t);
    /* * get current peak value */
    #[no_mangle]
    fn aubio_peakpicker_get_thresholded_input(p: *mut aubio_peakpicker_t)
     -> *mut fvec_t;
    /* * set peak picking threshold */
    #[no_mangle]
    fn aubio_peakpicker_set_threshold(p: *mut aubio_peakpicker_t,
                                      threshold: smpl_t) -> uint_t;
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
pub type smpl_t = f32;
/* * print format for sample in double precision */
/* * unsigned integer */
pub type uint_t = u32;
/* * signed integer */
pub type sint_t = i32;
/* * character */
pub type char_t = i8;
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
pub type C2RustUnnamed = u32;
pub const AUBIO_FAIL: C2RustUnnamed = 1;
pub const AUBIO_OK: C2RustUnnamed = 0;
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
/* * \file

  Logging features

  This file specifies ::aubio_log_set_function and
  ::aubio_log_set_level_function, which let you define one or several custom
  logging functions to redirect warnings and errors from aubio to your
  application. The custom function should have the prototype defined in
  ::aubio_log_function_t.

  After a call to ::aubio_log_set_level_function, ::aubio_log_reset can be used
  to reset each logging functions to the default ones.

  \example utils/test-log.c

*/
/* * list of logging levels */
pub type aubio_log_level = u32;
/* *< number of valid levels */
/* *< warnings */
pub const AUBIO_LOG_LAST_LEVEL: aubio_log_level = 5;
/* *< debug messages */
pub const AUBIO_LOG_WRN: aubio_log_level = 4;
/* *< general messages */
pub const AUBIO_LOG_DBG: aubio_log_level = 3;
/* *< infos */
pub const AUBIO_LOG_MSG: aubio_log_level = 2;
/* *< critical errors */
pub const AUBIO_LOG_INF: aubio_log_level = 1;
pub const AUBIO_LOG_ERR: aubio_log_level = 0;
#[derive(Copy, Clone)]
#[repr(C)]
pub struct cvec_t {
    pub length: uint_t,
    pub norm: *mut smpl_t,
    pub phas: *mut smpl_t,
}
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

  Spectral description functions

  All of the following spectral description functions take as arguments the FFT
  of a windowed signal (as created with aubio_pvoc). They output one smpl_t per
  buffer (stored in a vector of size [1]).

  \section specdesc Spectral description functions

  A list of the spectral description methods currently available follows.

  \subsection onsetdesc Onset detection functions

  These functions are designed to raise at notes attacks in music signals.

  \b \p energy : Energy based onset detection function

  This function calculates the local energy of the input spectral frame.

  \b \p hfc : High Frequency Content onset detection function

  This method computes the High Frequency Content (HFC) of the input spectral
  frame. The resulting function is efficient at detecting percussive onsets.

  Paul Masri. Computer modeling of Sound for Transformation and Synthesis of
  Musical Signal. PhD dissertation, University of Bristol, UK, 1996.

  \b \p complex : Complex Domain Method onset detection function

  Christopher Duxbury, Mike E. Davies, and Mark B. Sandler. Complex domain
  onset detection for musical signals. In Proceedings of the Digital Audio
  Effects Conference, DAFx-03, pages 90-93, London, UK, 2003.

  \b \p phase : Phase Based Method onset detection function

  Juan-Pablo Bello, Mike P. Davies, and Mark B. Sandler. Phase-based note onset
  detection for music signals. In Proceedings of the IEEE International
  Conference on Acoustics Speech and Signal Processing, pages 441­444,
  Hong-Kong, 2003.

  \b \p wphase : Weighted Phase Deviation onset detection function

  S. Dixon. Onset detection revisited. In Proceedings of the 9th International
  Conference on Digital Audio Ef- fects (DAFx) , pages 133–137, 2006.

  http://www.eecs.qmul.ac.uk/~simond/pub/2006/dafx.pdf

  \b \p specdiff : Spectral difference method onset detection function

  Jonhatan Foote and Shingo Uchihashi. The beat spectrum: a new approach to
  rhythm analysis. In IEEE International Conference on Multimedia and Expo
  (ICME 2001), pages 881­884, Tokyo, Japan, August 2001.

  \b \p kl : Kullback-Liebler onset detection function

  Stephen Hainsworth and Malcom Macleod. Onset detection in music audio
  signals. In Proceedings of the International Computer Music Conference
  (ICMC), Singapore, 2003.

  \b \p mkl : Modified Kullback-Liebler onset detection function

  Paul Brossier, ``Automatic annotation of musical audio for interactive
  systems'', Chapter 2, Temporal segmentation, PhD thesis, Centre for Digital
  music, Queen Mary University of London, London, UK, 2006.

  \b \p specflux : Spectral Flux

  Simon Dixon, Onset Detection Revisited, in ``Proceedings of the 9th
  International Conference on Digital Audio Effects'' (DAFx-06), Montreal,
  Canada, 2006.

  \subsection shapedesc Spectral shape descriptors

  The following descriptors are described in:

  Geoffroy Peeters, <i>A large set of audio features for sound description
  (similarity and classification) in the CUIDADO project</i>, CUIDADO I.S.T.
  Project Report 2004 (<a
  href="http://www.ircam.fr/anasyn/peeters/ARTICLES/Peeters_2003_cuidadoaudiofeatures.pdf">pdf</a>)

  \b \p centroid : Spectral centroid

  The spectral centroid represents the barycenter of the spectrum.

  \e Note: This function returns the result in bin. To get the spectral
  centroid in Hz, aubio_bintofreq() should be used.

  \b \p spread : Spectral spread

  The spectral spread is the variance of the spectral distribution around its
  centroid.

  See also <a href="http://en.wikipedia.org/wiki/Standard_deviation">Standard
  deviation</a> on Wikipedia.

  \b \p skewness : Spectral skewness

  Similarly, the skewness is computed from the third order moment of the
  spectrum. A negative skewness indicates more energy on the lower part of the
  spectrum. A positive skewness indicates more energy on the high frequency of
  the spectrum.

  See also <a href="http://en.wikipedia.org/wiki/Skewness">Skewness</a> on
  Wikipedia.

  \b \p kurtosis : Spectral kurtosis

  The kurtosis is a measure of the flatness of the spectrum, computed from the
  fourth order moment.

  See also <a href="http://en.wikipedia.org/wiki/Kurtosis">Kurtosis</a> on
  Wikipedia.

  \b \p slope : Spectral slope

  The spectral slope represents decreasing rate of the spectral amplitude,
  computed using a linear regression.

  \b \p decrease : Spectral decrease

  The spectral decrease is another representation of the decreasing rate,
  based on perceptual criteria.

  \b \p rolloff : Spectral roll-off

  This function returns the bin number below which 95% of the spectrum energy
  is found.

  \example spectral/test-specdesc.c

*/
/*
  Copyright (C) 2003-2015 Matthew Davies and Paul Brossier <piem@aubio.org>

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

  Beat tracking using a context dependant model

  This file implements the causal beat tracking algorithm designed by Matthew
  Davies and described in the following articles:

  Matthew E. P. Davies and Mark D. Plumbley. Causal tempo tracking of audio.
  In Proceedings of the International Symposium on Music Information Retrieval
  (ISMIR), pages 164­169, Barcelona, Spain, 2004.

  Matthew E. P. Davies, Paul Brossier, and Mark D. Plumbley. Beat tracking
  towards automatic musical accompaniment. In Proceedings of the Audio
  Engineering Society 118th Convention, Barcelona, Spain, May 2005.

  \example tempo/test-beattracking.c

*/
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

  Phase vocoder object

  This object implements a phase vocoder. The spectral frames are computed
  using a HanningZ window and a swapped version of the signal to simplify the
  phase relationships across frames. The window sizes and overlap are specified
  at creation time.

  \example spectral/test-phasevoc.c

*/
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

  Peak picking utilities function

  \example onset/test-peakpicker.c

*/
/*
  Copyright (C) 2006-2009 Paul Brossier <piem@aubio.org>

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
/* structure to store object state */
#[derive(Copy, Clone)]
#[repr(C)]
pub struct _aubio_tempo_t {
    pub od: *mut aubio_specdesc_t,
    pub pv: *mut aubio_pvoc_t,
    pub pp: *mut aubio_peakpicker_t,
    pub bt: *mut aubio_beattracking_t,
    pub fftgrain: *mut cvec_t,
    pub of: *mut fvec_t,
    pub dfframe: *mut fvec_t,
    pub out: *mut fvec_t,
    pub onset: *mut fvec_t,
    pub silence: smpl_t,
    pub threshold: smpl_t,
    pub blockpos: sint_t,
    pub winlen: uint_t,
    pub step: uint_t,
    pub samplerate: uint_t,
    pub hop_size: uint_t,
    pub total_frames: uint_t,
    pub last_beat: uint_t,
    pub delay: sint_t,
    pub last_tatum: uint_t,
    pub tatum_signature: uint_t,
}
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

  Tempo detection object

  This object stores all the memory required for tempo detection algorithm
  and returns the estimated beat locations.

  \example tempo/test-tempo.c
  \example examples/aubiotrack.c

*/
/* * tempo detection structure */
pub type aubio_tempo_t = _aubio_tempo_t;
/* * execute tempo detection

  \param o beat tracking object
  \param input new samples
  \param tempo output beats

*/
/* * number of tatum between each beats */
/* execute tempo detection function on iput buffer */
#[no_mangle]
pub unsafe extern "C" fn aubio_tempo_do(mut o: *mut aubio_tempo_t,
                                        mut input: *const fvec_t,
                                        mut tempo: *mut fvec_t) {
    let mut i: uint_t = 0;
    let mut winlen: uint_t = (*o).winlen;
    let mut step: uint_t = (*o).step;
    let mut thresholded: *mut fvec_t = 0 as *mut fvec_t;
    aubio_pvoc_do((*o).pv, input, (*o).fftgrain);
    aubio_specdesc_do((*o).od, (*o).fftgrain, (*o).of);
    /*if (usedoubled) {
    aubio_specdesc_do(o2,fftgrain, onset2);
    onset->data[0] *= onset2->data[0];
  }*/
  /* execute every overlap_size*step */
    if (*o).blockpos == step as i32 - 1 as i32 {
        /* check dfframe */
        aubio_beattracking_do((*o).bt, (*o).dfframe, (*o).out);
        /* rotate dfframe */
        i = 0 as i32 as uint_t;
        while i < winlen.wrapping_sub(step) {
            *(*(*o).dfframe).data.offset(i as isize) =
                *(*(*o).dfframe).data.offset(i.wrapping_add(step) as isize);
            i = i.wrapping_add(1)
        }
        i = winlen.wrapping_sub(step);
        while i < winlen {
            *(*(*o).dfframe).data.offset(i as isize) = 0.0f64 as smpl_t;
            i = i.wrapping_add(1)
        }
        (*o).blockpos = -(1 as i32)
    }
    (*o).blockpos += 1;
    aubio_peakpicker_do((*o).pp, (*o).of, (*o).onset);
    // store onset detection function in second sample of vector
  //tempo->data[1] = o->onset->data[0];
    thresholded = aubio_peakpicker_get_thresholded_input((*o).pp);
    *(*(*o).dfframe).data.offset(winlen.wrapping_sub(step).wrapping_add((*o).blockpos
                                                                            as
                                                                            u32)
                                     as isize) =
        *(*thresholded).data.offset(0 as i32 as isize);
    /* end of second level loop */
    *(*tempo).data.offset(0 as i32 as isize) =
        0 as i32 as smpl_t; /* reset tactus */
    //i=0;
    i = 1 as i32 as uint_t;
    while (i as f32) <
              *(*(*o).out).data.offset(0 as i32 as isize) {
        /* if current frame is a predicted tactus */
        if (*o).blockpos as f32 ==
               floorf(*(*(*o).out).data.offset(i as isize)) {
            *(*tempo).data.offset(0 as i32 as isize) =
                *(*(*o).out).data.offset(i as isize) -
                    floorf(*(*(*o).out).data.offset(i as
                                                        isize)); /* set tactus */
            /* test for silence */
            if aubio_silence_detection(input, (*o).silence) ==
                   1 as i32 as u32 {
                *(*tempo).data.offset(0 as i32 as isize) =
                    0 as i32 as smpl_t
                // unset beat if silent
            }
            (*o).last_beat =
                (*o).total_frames.wrapping_add(floorf(((*(*tempo).data.offset(0
                                                                                  as
                                                                                  i32
                                                                                  as
                                                                                  isize)
                                                            *
                                                            (*o).hop_size as
                                                                f32)
                                                           as f64 +
                                                           0.5f64) as
                                                          f32) as
                                                   uint_t);
            (*o).last_tatum = (*o).last_beat
        }
        i = i.wrapping_add(1)
    }
    (*o).total_frames =
        ((*o).total_frames as u32).wrapping_add((*o).hop_size) as
            uint_t as uint_t;
}
/* * get the time of the latest beat detected, in samples

  \param o tempo detection object as returned by ::new_aubio_tempo

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_tempo_get_last(mut o: *mut aubio_tempo_t)
 -> uint_t {
    return (*o).last_beat.wrapping_add((*o).delay as u32);
}
/* * get the time of the latest beat detected, in seconds

  \param o tempo detection object as returned by ::new_aubio_tempo

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_tempo_get_last_s(mut o: *mut aubio_tempo_t)
 -> smpl_t {
    return aubio_tempo_get_last(o) as f32 /
               (*o).samplerate as smpl_t;
}
/* * get the time of the latest beat detected, in milliseconds

  \param o tempo detection object as returned by ::new_aubio_tempo

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_tempo_get_last_ms(mut o: *mut aubio_tempo_t)
 -> smpl_t {
    return (aubio_tempo_get_last_s(o) as f64 * 1000.0f64) as
               smpl_t;
}
/* * set current delay

  \param o beat tracking object
  \param delay delay to set tempo to, in samples

  \return `0` if successful, non-zero otherwise

 */
#[no_mangle]
pub unsafe extern "C" fn aubio_tempo_set_delay(mut o: *mut aubio_tempo_t,
                                               mut delay: sint_t) -> uint_t {
    (*o).delay = delay;
    return AUBIO_OK as i32 as uint_t;
}
/* * set current delay in seconds

  \param o beat tracking object
  \param delay delay to set tempo to, in seconds

  \return `0` if successful, non-zero otherwise

 */
#[no_mangle]
pub unsafe extern "C" fn aubio_tempo_set_delay_s(mut o: *mut aubio_tempo_t,
                                                 mut delay: smpl_t)
 -> uint_t {
    (*o).delay = (delay * (*o).samplerate as f32) as sint_t;
    return AUBIO_OK as i32 as uint_t;
}
/* * set current delay

  \param o beat tracking object
  \param delay delay to set tempo to, in samples

  \return `0` if successful, non-zero otherwise

 */
#[no_mangle]
pub unsafe extern "C" fn aubio_tempo_set_delay_ms(mut o: *mut aubio_tempo_t,
                                                  mut delay: smpl_t)
 -> uint_t {
    return aubio_tempo_set_delay_s(o,
                                   (delay as f64 / 1000.0f64) as
                                       smpl_t);
}
/* * get current delay

  \param o beat tracking object

  \return current delay, in samples

 */
#[no_mangle]
pub unsafe extern "C" fn aubio_tempo_get_delay(mut o: *mut aubio_tempo_t)
 -> uint_t {
    return (*o).delay as uint_t;
}
/* * get current delay in seconds

  \param o beat tracking object

  \return current delay, in seconds

 */
#[no_mangle]
pub unsafe extern "C" fn aubio_tempo_get_delay_s(mut o: *mut aubio_tempo_t)
 -> smpl_t {
    return (*o).delay as f32 / (*o).samplerate as smpl_t;
}
/* * get current delay in ms

  \param o beat tracking object

  \return current delay, in milliseconds

 */
#[no_mangle]
pub unsafe extern "C" fn aubio_tempo_get_delay_ms(mut o: *mut aubio_tempo_t)
 -> smpl_t {
    return (aubio_tempo_get_delay_s(o) as f64 * 1000.0f64) as
               smpl_t;
}
/* * set tempo detection silence threshold

  \param o beat tracking object
  \param silence new silence threshold, in dB

  \return `0` if successful, non-zero otherwise

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_tempo_set_silence(mut o: *mut aubio_tempo_t,
                                                 mut silence: smpl_t)
 -> uint_t {
    (*o).silence = silence;
    return AUBIO_OK as i32 as uint_t;
}
/* * get tempo detection silence threshold

  \param o tempo detection object as returned by new_aubio_tempo()

  \return current silence threshold

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_tempo_get_silence(mut o: *mut aubio_tempo_t)
 -> smpl_t {
    return (*o).silence;
}
/* * set tempo detection peak picking threshold

  \param o beat tracking object
  \param threshold new threshold

  \return `0` if successful, non-zero otherwise

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_tempo_set_threshold(mut o: *mut aubio_tempo_t,
                                                   mut threshold: smpl_t)
 -> uint_t {
    (*o).threshold = threshold;
    aubio_peakpicker_set_threshold((*o).pp, (*o).threshold);
    return AUBIO_OK as i32 as uint_t;
}
/* * get tempo peak picking threshold

  \param o tempo detection object as returned by new_aubio_tempo()

  \return current tempo detection threshold

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_tempo_get_threshold(mut o: *mut aubio_tempo_t)
 -> smpl_t {
    return (*o).threshold;
}
/* * create tempo detection object

  \param method beat tracking method, unused for now (use "default")
  \param buf_size length of FFT
  \param hop_size number of frames between two consecutive runs
  \param samplerate sampling rate of the signal to analyze

  \return newly created ::aubio_tempo_t if successful, `NULL` otherwise

*/
/* Allocate memory for an tempo detection */
#[no_mangle]
pub unsafe extern "C" fn new_aubio_tempo(mut tempo_mode: *const char_t,
                                         mut buf_size: uint_t,
                                         mut hop_size: uint_t,
                                         mut samplerate: uint_t)
 -> *mut aubio_tempo_t {
    let mut o: *mut aubio_tempo_t =
        calloc(::std::mem::size_of::<aubio_tempo_t>() as u64,
               1 as i32 as u64) as *mut aubio_tempo_t;
    let mut specdesc_func: [char_t; 1024] = [0; 1024];
    (*o).samplerate = samplerate;
    // check parameters are valid
    if (hop_size as sint_t) < 1 as i32 {
        aubio_log(AUBIO_LOG_ERR as i32,
                  b"AUBIO ERROR: tempo: got hop size %d, but can not be < 1\n\x00"
                      as *const u8 as *const i8, hop_size);
    } else if (buf_size as sint_t) < 2 as i32 {
        aubio_log(AUBIO_LOG_ERR as i32,
                  b"AUBIO ERROR: tempo: got window size %d, but can not be < 2\n\x00"
                      as *const u8 as *const i8, buf_size);
    } else if buf_size < hop_size {
        aubio_log(AUBIO_LOG_ERR as i32,
                  b"AUBIO ERROR: tempo: hop size (%d) is larger than window size (%d)\n\x00"
                      as *const u8 as *const i8, buf_size,
                  hop_size);
    } else if (samplerate as sint_t) < 1 as i32 {
        aubio_log(AUBIO_LOG_ERR as i32,
                  b"AUBIO ERROR: tempo: samplerate (%d) can not be < 1\n\x00"
                      as *const u8 as *const i8, samplerate);
    } else {
        /* length of observations, worth about 6 seconds */
        (*o).winlen =
            aubio_next_power_of_two((5.8f64 * samplerate as f64 /
                                         hop_size as f64) as
                                        uint_t);
        if (*o).winlen < 4 as i32 as u32 {
            (*o).winlen = 4 as i32 as uint_t
        }
        (*o).step =
            (*o).winlen.wrapping_div(4 as i32 as u32);
        (*o).blockpos = 0 as i32;
        (*o).threshold = 0.3f64 as smpl_t;
        (*o).silence = -90.0f64 as smpl_t;
        (*o).total_frames = 0 as i32 as uint_t;
        (*o).last_beat = 0 as i32 as uint_t;
        (*o).delay = 0 as i32;
        (*o).hop_size = hop_size;
        (*o).dfframe = new_fvec((*o).winlen);
        (*o).fftgrain = new_cvec(buf_size);
        (*o).out = new_fvec((*o).step);
        (*o).pv = new_aubio_pvoc(buf_size, hop_size);
        (*o).pp = new_aubio_peakpicker();
        aubio_peakpicker_set_threshold((*o).pp, (*o).threshold);
        if strcmp(tempo_mode,
                  b"default\x00" as *const u8 as *const i8) ==
               0 as i32 {
            strncpy(specdesc_func.as_mut_ptr(),
                    b"specflux\x00" as *const u8 as *const i8,
                    (1024 as i32 - 1 as i32) as
                        u64);
        } else {
            strncpy(specdesc_func.as_mut_ptr(), tempo_mode,
                    (1024 as i32 - 1 as i32) as
                        u64);
            specdesc_func[(1024 as i32 - 1 as i32) as usize] =
                '\u{0}' as i32 as char_t
        }
        (*o).od = new_aubio_specdesc(specdesc_func.as_mut_ptr(), buf_size);
        (*o).of = new_fvec(1 as i32 as uint_t);
        (*o).bt =
            new_aubio_beattracking((*o).winlen, (*o).hop_size,
                                   (*o).samplerate);
        (*o).onset = new_fvec(1 as i32 as uint_t);
        /*if (usedoubled)    {
    o2 = new_aubio_specdesc(type_onset2,buffer_size);
    onset2 = new_fvec(1);
  }*/
        if (*o).dfframe.is_null() || (*o).fftgrain.is_null() ||
               (*o).out.is_null() || (*o).pv.is_null() || (*o).pp.is_null() ||
               (*o).od.is_null() || (*o).of.is_null() || (*o).bt.is_null() ||
               (*o).onset.is_null() {
            aubio_log(AUBIO_LOG_ERR as i32,
                      b"AUBIO ERROR: tempo: failed creating tempo object\n\x00"
                          as *const u8 as *const i8);
        } else {
            (*o).last_tatum = 0 as i32 as uint_t;
            (*o).tatum_signature = 4 as i32 as uint_t;
            return o
        }
    }
    del_aubio_tempo(o);
    return 0 as *mut aubio_tempo_t;
}
/* * get current tempo

  \param o beat tracking object

  \return the currently observed tempo, or `0` if no consistent value is found

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_tempo_get_bpm(mut o: *mut aubio_tempo_t)
 -> smpl_t {
    return aubio_beattracking_get_bpm((*o).bt);
}
/* * get current beat period in samples

  \param bt beat tracking object

  Returns the currently observed period, in samples, or 0 if no consistent
  value is found.

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_tempo_get_period(mut o: *mut aubio_tempo_t)
 -> smpl_t {
    return aubio_beattracking_get_period((*o).bt);
}
/* * get current beat period in seconds

  \param bt beat tracking object

  Returns the currently observed period, in seconds, or 0 if no consistent
  value is found.

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_tempo_get_period_s(mut o: *mut aubio_tempo_t)
 -> smpl_t {
    return aubio_beattracking_get_period_s((*o).bt);
}
/* * get current tempo confidence

  \param o beat tracking object

  \return confidence with which the tempo has been observed, the higher the
  more confidence, `0` if no consistent value is found.

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_tempo_get_confidence(mut o: *mut aubio_tempo_t)
 -> smpl_t {
    return aubio_beattracking_get_confidence((*o).bt);
}
/* * check whether a tatum was detected in the current frame

   \param o beat tracking object

   \return 2 if a beat was detected, 1 if a tatum was detected, 0 otherwise

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_tempo_was_tatum(mut o: *mut aubio_tempo_t)
 -> uint_t {
    let mut last_tatum_distance: uint_t =
        (*o).total_frames.wrapping_sub((*o).last_tatum);
    let mut beat_period: smpl_t = aubio_tempo_get_period(o);
    let mut tatum_period: smpl_t =
        beat_period / (*o).tatum_signature as f32;
    if last_tatum_distance < (*o).hop_size {
        (*o).last_tatum = (*o).last_beat;
        return 2 as i32 as uint_t
    } else {
        if last_tatum_distance as f32 > tatum_period {
            if last_tatum_distance.wrapping_add((*o).hop_size) as
                   f32 > beat_period {
                // next beat is too close, pass
                return 0 as i32 as uint_t
            }
            (*o).last_tatum = (*o).total_frames;
            return 1 as i32 as uint_t
        }
    }
    return 0 as i32 as uint_t;
}
/* * get position of last_tatum, in samples

   \param o beat tracking object

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_tempo_get_last_tatum(mut o: *mut aubio_tempo_t)
 -> smpl_t {
    return (*o).last_tatum as smpl_t - (*o).delay as f32;
}
/* * set number of tatum per beat

   \param o beat tracking object
   \param signature number of tatum per beat (between 1 and 64)

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_tempo_set_tatum_signature(mut o:
                                                             *mut aubio_tempo_t,
                                                         mut signature:
                                                             uint_t)
 -> uint_t {
    if signature < 1 as i32 as u32 ||
           signature > 64 as i32 as u32 {
        return AUBIO_FAIL as i32 as uint_t
    } else {
        (*o).tatum_signature = signature;
        return AUBIO_OK as i32 as uint_t
    };
}
/* * delete tempo detection object

  \param o beat tracking object

*/
#[no_mangle]
pub unsafe extern "C" fn del_aubio_tempo(mut o: *mut aubio_tempo_t) {
    if !(*o).od.is_null() { del_aubio_specdesc((*o).od); }
    if !(*o).bt.is_null() { del_aubio_beattracking((*o).bt); }
    if !(*o).pp.is_null() { del_aubio_peakpicker((*o).pp); }
    if !(*o).pv.is_null() { del_aubio_pvoc((*o).pv); }
    if !(*o).out.is_null() { del_fvec((*o).out); }
    if !(*o).of.is_null() { del_fvec((*o).of); }
    if !(*o).fftgrain.is_null() { del_cvec((*o).fftgrain); }
    if !(*o).dfframe.is_null() { del_fvec((*o).dfframe); }
    if !(*o).onset.is_null() { del_fvec((*o).onset); }
    free(o as *mut core::ffi::c_void);
}
