;fuction that returns (abs(z))^2

function modsq, z

return, real_part(z*conj(z))
end
