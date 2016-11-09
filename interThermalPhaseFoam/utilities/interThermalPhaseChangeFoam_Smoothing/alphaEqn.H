scalar smoothingFactor = 0.5;
//pass1
alpha1 = smoothingFactor*(fvc::average(fvc::interpolate(alpha1))) + (1.0 - smoothingFactor)*alpha1;
//pass2
alpha1 = smoothingFactor*(fvc::average(fvc::interpolate(alpha1))) + (1.0 - smoothingFactor)*alpha1;
//pass3
alpha1 = smoothingFactor*(fvc::average(fvc::interpolate(alpha1))) + (1.0 - smoothingFactor)*alpha1;

