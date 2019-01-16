//computes cumulative length of r1 and r2;
size_t N = r1.size() + r2.size() - 1;
size_t M = fft_next_good_size(N);

//padds signals with zeros
while (r1.size() != M)
    r1.emplace_back(0.0);
while (r2.size() != M)
    r2.emplace_back(0.0);

//Performs FWD FFT1
fft.fwd(fv1, r1);

//Performs FWD FFT2
fft.fwd(fv2, r2);

//multiply FV1 * FV2, element-wise
std::transform(fv1.begin(), fv1.end(),
               fv2.begin(), std::back_inserter(mulvec),
               std::multiplies<complex<float>>());

//FFT INVERSE
fft.inv(result, mulvec);
return result;