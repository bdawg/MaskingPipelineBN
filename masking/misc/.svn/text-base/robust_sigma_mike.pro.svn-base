;;Robust sigma repeatedly gives dumb results, e.g. for a 2-valued distribution.
function robust_sigma_mike, y
ny = n_elements(y)
medy = median(y)
;;Find the squared deviations from the median rather than the mean
;;(slightly more robust)
square_devs = (y-medy)^2
s=sort(square_devs)
;;1 sigma is the 68th percentile for a Gaussian distribution. This is
;;really not "correct" for a small number of inputs. But in that case, why 
;;are you calculating sigma, as you have a student's t-distribution with
;;giant wings? i.e. this is just as dodgy as assuming a Gaussian error.
return, sqrt(square_devs[s[0.68*ny]])
end
