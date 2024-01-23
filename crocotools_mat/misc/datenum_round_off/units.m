function x = units( x, unit1, unit2 )
%
% y = units(x,'unit1','unit2')  - Converts x from physical unit1 to unit2
%
% example: units(1,'in','mm') = 25.4
%
% The following units are supported:
% Acceleration: m/s^2, cm/s^2, mm/s^2, ft/s^2, in/s^2, G
% Angle: rad, deg, rev
% Area: km^2, m^2, cm^2, mm^2, ym^2 (square-micrometer), sq-mile, sq-yd, sq-ft, sq-in,
%    acres, ha, ar
% Area Moment of Inertia: m^4, cm^4, mm^4, ft^4, in^4
% Density: t/m^3 (metric), kg/m^3, g/cm^3, g/mm^3, lbs/ft^3 lbs/in^3, lbs/galUS, lbs/galUK
% Energy, Work, Torque: GJ, MJ, kJ, J, mJ, Nm, Ncm, Nmm, kWh, Wh, Ws, lb-ft, lb-in, oz-in,Btu,
%    Btu, cal, kcal, eV
% Force: MN, kN, N, dyne, lbf, kip
% Fuel Consumption: l/100km, miles/galUS
% Frequency, Angular Velocity: GHz, MHz, kHz, Hz, 1/s, 1/min, 1/h, rad/s, deg/s, rpm
% Length: km, m, dm, cm, mm, ym (micrometer), nm, mile, yard, ft, in, mill,
%    Angstrom, light-year
% Mass: t (metric), tUS, tUK, kg, g, mg, yg (microgram), ng, lbs, oz
% Mass Moment of Inertia: kg*m^2, kg*cm^2, kg*mm^2, g*m^2, g*cm^2, g*mm^2, lb*ft^2, lb*in^2
% Power: GW, MW, kW, W, mW, hp, Btu/h, Btu/s, kcal/h, J/h
% Pressure, Stress: GPa, MPa, kPa, hPa, Pa, bar, mbar, atm, dyne/cm^2, ksi, psi, mmHg, mmH2O
% Strain: m/m, mm/m, ym/m (micrometer/m), nm/m, %, o/oo, in/in, mill/in
% Stress Intensity Factor: MPa*m^1/2, MPa*mm^1/2, ksi*in^1/2, psi*in^1/2
% Temperature: degK, degC, degF, degR
% Time: yr (365 days), mth (30 days), wk, day, hr, min, s, ms, ys (microsecond), ns
% Velocity: km/s, km/h, m/s, cm/s, mm/s, m/min, mm/min, mps, mph, ft/s, in/s, ft/min,
%    in/min, Mach, knots
% Viscosity: Ns/m^2, poise, centipoise, lbfs/sq-ft
% Volume: km^3, m^3, cm^3, mm^3, ym^3, cu-mile, cu-ft, cu-in, l, cl, ml, galUS, galUK,
%    pint (liquid US), quart (liquid US), fl-oz (liquid US)
%
% Metric Prefixes: Yotta, Zetta, Exa, Peta, Terra, Giga, Mega, Myria, kilo, hecto, 1
%     deci, centi, milli, micro, nano, pico, femto, atto, zepto, yocto
%
%
% CONSISTENCY OF UNITS IS NOT CHECKED FOR!
% (i.e., units(1,'m','kg') = 1, which is of course nonsense.)
%


% -------------------------------------------------------------------------
% Henning Ressing, PhD
% The University of British Columbia, Mechanical Engineering
% Last updated: November 9, 2003



err = 0;

switch unit1

% Length
    case 'km'
      x = x.*1000;  
    case {'m', 'meter', 'meters'}
      %x = x;
    case 'dm'
      x = x./10;
    case 'cm'
      x = x./100;
    case 'mm'
      x = x./1000;
    case 'ym'
      x = x./1e6;
    case 'nm'
      x = x./1e9;
    case {'mile', 'miles'}
      x = x.*1609.344;
    case {'nautical mile', 'nautical miles'}
      x = x.*1852;
    case {'yd', 'yard', 'yards'}
      x = x.*0.9144;
    case {'f', 'ft', 'feet', 'foot'}
      x = x.*0.3048;  
    case {'in', 'inch', 'inches'}
      x = x.*0.0254;
    case 'mill'
      x = x.*(25.4e-6);
    case 'Angstrom'
      x = x./1e10;
    case {'light-year', 'light-yr'}
      x = x.*9.46075309081898e15;

    % Time
    case {'yr', 'year', 'years'}
      x = x.*(60*60*24*365);  
    case {'mth', 'mon', 'month', 'months'}
      x = x.*(60*60*24*30);
    case {'wk', 'week', 'weeks'}
      x = x.*(60*60*24*7);
    case {'day', 'days'}
      x = x.*(60*60*24);
    case {'hr', 'hour', 'hours'}
      x = x.*(60*60);
    case {'min', 'minute', 'minutes'}
      x = x.*60;
    case {'s', 'sec', 'second', 'seconds'}
      %x = x;
    case {'ms', 'msec', 'millisec', 'millisecond', 'milliseconds'}
      x = x./1000;
    case 'ys'
      x = x./1e6;
    case {'ns', 'nanosec', 'nanosecond', 'nanoseconds'}
      x = x./1e9;

    % Mass
    case 't'
      x = x.*1000;
    case 'tUS'
      x = x.*907.18474;
    case 'tUK'
      x = x.*1016.0469088;
    case {'kg', 'kilogram', 'kilograms'}
      %x = x;  
    case {'g', 'gram', 'gm', 'grams'}
      x = x./1000;
    case {'mg', 'milligram', 'milligrams'}
      x = x./1e6;
    case 'yg'
      x = x./1e9;
    case 'ng'
      x = x./1e12;
    case 'lbs'
      x = x.*0.4535924;
    case {'oz', 'ounce', 'ounces'}
      x = x.*0.03110348;

    % Force
    case 'MN'
      x = x.*1e6;  
    case 'kN'
      x = x.*1000;  
    case {'N', 'Newton', 'newton', 'Newtons', 'newtons'}
      %x = x;
    case 'dyne'
      x = x./1e5;
    case 'lbf'
      x = x.*4.448222;
    case 'kip'
      x = x.*4448.222;

    % Pressure, Stress
    case 'GPa'
      x = x.*1e9;  
    case 'MPa'
      x = x.*1e6;
    case 'N/mm^2'
      x = x.*1e6;
    case 'kPa'
      x = x.*1000;
    case 'hPa'
      x = x.*100;
    case 'Pa'
      %x = x;
    case 'bar'
      x = x.*1e5;
    case 'mbar'
      x = x.*100;
    case 'dyne/cm^2'
      x = x./10;
    case 'atm'
      x = x.*1.01325e5;
    case 'ksi'
      x = x.*6894757;
    case 'psi'
      x = x.*6894.757;
    case 'mmHg'
      x = x.*133.322;
    case 'mmH2O'
      x = x.*9.80665;

    % Temperature
    case 'degK'
      %x = x;
    case 'degC'
      x = x + 273.15;
    case 'degF'
      x = 5/9.*(x-32) + 273.15;
    case 'degR'
      x = 5/9.*(x-491.67) + 273.15;

    % Work, Energy, Torque
    case 'GJ'
      x = x.*10^9;
    case 'MJ'
      x = x.*10^6;
    case 'kJ'
      x = x.*1000;
    case 'J'
      %x = x;
    case 'mJ'
      x = x./1000;
    case 'Nm'
      %x = x;
    case 'Ncm'
      x = x./100;
    case 'Nmm'
      x = x./1000;
    case 'lb-ft'
      x = x.*1.355818;
    case 'lb-in'
      x = x.*0.112984833;
    case 'oz-in'
      x = x.*7.061552e-3;
    case 'Btu'
      x = x.*1055.056;
    case 'kWh'
      x = x.*3.6e6;
    case 'Wh'
      x = x.*3600;
    case 'Ws'
      %x = x;
    case 'cal'
      x = x.*4.1868;
    case 'kcal'
      x = x.*4186.8;
    case 'eV'
      x = x.*1.60218e-19;

    % Power
    case 'GW'
      x = x.*1e9;
    case 'MW'
      x = x.*1e6;
    case 'kW'
      x = x.*1e3;
    case 'W'
      %x = x;
    case 'mW'
      x = x./1000;
    case 'hp'
      x = x.*745.6999;
    case 'Btu/h'
      x = x.*0.2930711;
    case 'Btu/s'
      x = x.*1055.056;
    case 'kcal/h'
      x = x.*1.163;
    case 'J/h'
      x = x./3600;

    % Velocity
    case {'km/h', 'km/hr'}
      x = x.*(1000/3600);
    case {'km/s', 'km/sec'}
      x = x.*1000;
    case {'m/s', 'm/sec', 'meter/sec', 'meters/sec'}
      %x = x;
    case {'cm/s', 'cm/sec'}
      x = x./100;
    case 'mm/s'
      x = x./1000;
    case 'm/min'
      x = x./60;
    case 'mm/min'
      x = x./60000;
    case 'mph'
      x = x.*(1609.344/3600);
    case 'mps'
      x = x.*1609.344;
    case {'ft/s', 'ft/sec'}
      x = x.*0.3048;
    case {'in/s', 'in/sec'}
      x = x.*0.0254;
    case 'ft/min'
      x = x.*(0.3048/60);
    case 'in/min'
      x = x.*(0.0254/60);
    case 'Mach'
      x = x.*331.5;
    case {'kts', 'knot', 'knots'}
      x = x.*0.514444;

    % Acceleration
    case 'm/s^2'
      %x = x;
    case 'cm/s^2'
      x = x./100;
    case 'mm/s^2'
      x = x./1000;
    case 'ft/s^2'
      x = x.*0.3048;
    case 'in/s^2'
      x = x.*0.0254;
    case 'G'
      x = x.*9.81;

    % Area
    case 'km^2'
      x = x.*1e6;
    case 'm^2'
      %x = x;
    case 'cm^2'
      x = x./1e4;
    case 'mm^2'
      x = x./1e6;
    case 'ym^2'
      x = x./1e12;
    case 'sq-mile'
      x = x.*1609.344^2;
    case 'sq-yd'
      x = x.*0.83612736;
    case 'sq-ft'
      x = x.*0.3048^2;
    case 'sq-in'
      x = x.*0.0254^2;
    case 'acres'
      x = x.*4046.8564224;
    case 'ha'
      x = x.*1e4;
    case 'ar'
      x = x.*100;

    % Volume
    case 'km^3'
      x = x.*10^9;
    case 'm^3'
      %x = x;
    case 'cm^3'
      x = x./10^6;
    case 'mm^3'
      x = x./10^9;
    case 'ym^3'
      x = x./10^18;
    case 'cu-mile'
      x = x.*1609.344^3;
    case 'cu-yd'
      x = x.*0.764554857984;
    case 'cu-ft'
      x = x.*0.3048^3;
    case 'cu-in'
      x = x.*0.0254^3;
    case 'l'
      x = x./1000;
    case 'cl'
      x = x./1e6;
    case 'ml'
      x = x./10^9;
    case 'galUS'
      x = x.*3.785412e-3;
    case 'galUK'
      x = x.*4.54609e-3; 
    case 'pint'
      x = x.*4.73176473e-4;
    case 'quart'
      x = x.*9.46352946e-4;
    case 'fl-oz'
      x = x.*2.95735295625e-5;

    % Density
    case 't//m^3'
      x = x.*1000;
    case 'kg/m^3'
      %x = x;
    case 'g/cm^3'
      x = x.*1000;
    case 'g/mm^3'
      x = x.*1e6;
    case 'kg/l'
      x = x.*1000;
    case 'lbs/ft^3'
      x = x.*(0.4535924/0.3048^3);
    case 'lbs/in^3'
      x = x.*(0.4535924/0.0254^3);
    case 'lbs/galUS'
      x = x.*119.826427;
    case 'lbs/galUK'
      x = x.*99.776373;

    % Mass Moment of Inertia
    case 'kg*m^2'
      %x = x;
    case 'kg*cm^2'
      x = x./1e4;
    case 'kg*mm^2'
      x = x./1e6;
    case 'g*m^2'
      x = x./1000;
    case 'g*cm^2'
      x = x./1e7;
    case 'g*mm^2'
      x = x./1e9;
    case 'lb*ft^2'
      x = x.*0.04214011;
    case 'lb*in^2'
      x = x.*0.2926397e-3;

    % Area Moment of Inertia
    case 'm^4'
      %x = x;
    case 'cm^4'
      x = x./1e8;
    case 'mm^4'
      x = x./1e12;
    case 'ft^4'
      x = x./(0.3048^4);
    case 'in^4'
      x = x./(0.0254^4);

    % Frequency / Angular Velocity
    case 'GHz'
      x = x.*1e9;
    case 'MHz'
      x = x.*1e6;
    case 'kHz'
      x = x.*1e3;
    case 'Hz'
      %x = x;
    case '1/min'
      x = x./60;
    case '1/h'
      x = x./3600;
    case 'rad/s'
      x = x./(2*pi);
    case 'deg/s'
      x = x./360;
    case 'rpm'
      x = x./60;

    % Angle
    case 'rad'
      %x = x;
    case 'deg'
      x = x.*(pi/180);
    case 'rev'
      x = x.*(2*pi);

    % Stress Intensity Factor
    case 'MPa*m^1/2'
      %x = x;
    case 'MPa*mm^1/2'
      x = x./sqrt(1000);
    case 'ksi*in^1/2'
      x = x.*6.894757.*sqrt(0.0254);
    case 'psi*in^1/2'
      x = x.*6894.757.*sqrt(0.0254);

    % Fuel Consumption
    case 'l/100km'
      %x = x;
    case 'miles/galUS'
      x = 235.214596754951./x;

    % Viscosity
    case 'Ns/m^2'
      %x = x;
    case 'poise'
      x = x./10;
    case 'centipoise'
      x = x./1000;
    case 'lbfs/sq-ft'
      x = x./0.02089;

    % Strain
    case 'm/m'
      %x = x;
    case 'mm/m'
      x = x./1e3;
    case 'ym/m'
      x = x./1e6;
    case 'nm/m'
      x = x./1e9;
    case '%'
      x = x./1e2;
    case 'o/oo'
      x = x./1e3;
    case 'in/in'
      %x = x;
    case 'mill/in'
      x = x./1e3';

    % Metric Prefixes
    case 'Yotta'
      x = x.*1e24;
    case 'Zetta'
      x = x.*1e21;
    case 'Exa'
      x = x.*1e18;
    case 'Peta'
      x = x.*1e15;
    case 'Tera'
      x = x.*1e12;
    case 'Giga'
      x = x.*1e9;
    case 'Mega'
      x = x.*1e6;
    case 'Myria'
      x = x.*1e5;
    case 'kilo'
      x = x.*1e3;
    case 'hecto'
      x = x.*1e2;
    case '1'
      %x = x;
    case 'deci'
      x = x.*1e-1;
    case 'centi'
      x = x.*1e-2;
    case 'milli'
      x = x.*1e-3;
    case 'micro'
      x = x.*1e-6;
    case 'nano'
      x = x.*1e-9;
    case 'pico'
      x = x.*1e-12;
    case 'femto'
      x = x.*1e-15;
    case 'atto'
      x = x.*1e-18;
    case 'zepto'
      x = x.*1e-21;
    case 'yocto'
      x = x.*1e-24;


otherwise
  err = 1;
  write_error_to_log(['Error - Unsupported unit1:', unit1], mfilename);

end

% ----------------------------------------------------------

if ~err

switch unit2
  
    % Length
    case {'km', 'kilometer', 'kilometers'}
      x = x./1000;  
    case {'m', 'meter', 'meters'}
      %x = x;
    case {'dm', 'decimeter', 'decimeters'}
      x = x.*10;
    case {'cm', 'centimeter', 'centimeters'}
      x = x.*100;
    case {'mm', 'millimeter', 'millimeters'}
      x = x.*1000;
    case 'ym'
      x = x.*1e6;
    case 'nm'
      x = x.*1e9;
    case {'mile', 'miles'}
      x = x./1609.344;
	case {'nautical mile', 'nautical miles'}
      x = x./1852;
    case {'yard', 'yd', 'yards'}
      x = x./0.9144;
    case {'f', 'ft', 'feet', 'foot'}
      x = x./0.3048;  
    case {'in', 'inch', 'inches'}
      x = x./0.0254;
    case 'mill'
      x = x./25.4e-6;
    case 'Angstrom'
      x = x.*1e10;
    case 'light-year'
      x = x./9.46075309081898e15;

    % Time
    case {'yr', 'year', 'years'}
      x = x./(60*60*24*365);  
    case {'mth', 'mon', 'month', 'months'}
      x = x./(60*60*24*30);
    case {'wk', 'week', 'weeks'}
      x = x./(60*60*24*7);
    case {'day', 'days'}
      x = x./(60*60*24);
    case {'hr', 'hour', 'hours'}
      x = x./(60*60);
    case {'min', 'minute', 'minutes'}
      x = x./60;
    case {'s', 'sec', 'second', 'seconds'}
      %x = x;
    case {'ms', 'millisec', 'millisecond', 'milliseconds'}
      x = x.*1000;
    case 'ys';
      x = x.*1e6;
    case 'ns';
      x = x.*1e9;

    % Mass
    case 't'
      x = x./1000;
    case 'tUS'
      x = x./907.18474;
    case 'tUK'
      x = x./1016.0469088;
    case {'kg', 'kilogram', 'kilograms'}
      %x = x;  
    case {'g', 'gram', 'grams'}
      x = x.*1000;
    case 'mg'
      x = x.*1e6;
    case 'yg'
      x = x.*1e9;
    case 'ng'
      x = x.*1e12;
    case 'lbs'
      x = x./0.4535924;
    case {'oz', 'ounce', 'ounces'} 
      x = x./0.03110348;

    % Force
    case 'MN'
      x = x./1e6;  
    case 'kN'
      x = x./1000;  
    case 'N'
      %%x = x;
    case 'dyne'
      x = x.*1e5;
    case 'lbf'
      x = x./4.44822;
    case 'kip'
      x = x./4448.222;

    % Pressure, Stress
    case 'GPa'
      x = x./10^9;  
    case 'MPa'
      x = x./10^6;  
    case 'N/mm^2'
      x = x./10^6;  
    case 'kPa'
      x = x./1000;
    case 'hPa'
      x = x./100;  
    case 'Pa'
      %%x = x;
    case 'bar'
      x = x./1e5;
    case 'mbar'
      x = x./100;
    case 'dyne/cm^2'
      x = x.*10;
    case 'atm'
      x = x./1.01325e5;
    case 'ksi'
      x = x./6894757;
    case 'psi'
      x = x./6894.757;
    case 'mmHg'
      x = x./133.322;
    case 'mmH2O'
      x = x./9.80665;

    % Temperature
    case 'degK'
      %%x = x;
    case 'degC'
      x = x - 273.15;
    case 'degF'
      x = 9/5.*(x-273.15) + 32;
    case 'degR'
      x = 9/5.*(x-273.15) + 491.67;

    % Work, Energy, Torque
    case 'GJ'
      x = x./10^9;
    case 'MJ'
      x = x./10^6;
    case 'kJ'
      x = x./1000;
    case 'J'
      %%x = x;
    case 'mJ'
      x = x.*1000;
    case 'Nm'
      %%x = x;
    case 'Ncm'
      x = x.*100;
    case 'Nmm'
      x = x.*1000;
    case 'lb-ft'
      x = x./1.355818;
    case 'lb-in'
      x = x./0.112984833;
    case 'oz-in'
      x = x./0.007061552;
    case 'Btu'
      x = x./1055.056;
    case 'kWh'
      x = x./3.6e6;
    case 'Wh'
      x = x./3600;
    case 'Ws'
      %x = x;
    case 'cal'
      x = x./4.1868;
    case 'kcal'
      x = x./4186.8;
    case 'eV'
      x = x./1.60218e-19;

    % Power
    case 'GW'
      x = x./1e9;
    case 'MW'
      x = x./1e6;
    case 'kW'
      x = x./1000;
    case 'W'
      %x = x;
    case 'mW'
      x = x.*1000;
    case 'hp'
      x = x./745.6999;
    case 'Btu/h'
      x = x./0.2930711;
    case 'Btu/s'
      x = x./1055.056;
    case 'kcal/h'
      x = x./1.163;
    case 'J/h'
      x = x.*3600;

    % Velocity
    case {'km/h', 'km/hr'}
      x = x./(1000/3600);
    case 'km/s'
      x = x./1000;
    case {'m/s', 'm/sec', 'meter/sec', 'meters/sec'}
      %x = x;
    case 'cm/s'
      x = x.*100;
    case 'mm/s'
      x = x.*1000;
    case 'm/min'
      x = x.*60;
    case 'mm/min'
      x = x.*60000;
    case 'mph'
      x = x./(1609.344/3600);
    case 'mps'
      x = x./1609.344;
    case {'ft/s', 'ft/sec'}
      x = x./0.3048;
    case 'in/s'
      x = x./0.0254;
    case 'ft/min'
      x = x./(0.3048*60);
    case 'in/min'
      x = x./(0.0254*60);
    case 'Mach'
      x = x./331.5;
    case {'kts', 'knot', 'knots'}
      x = x./0.514444;

    % Acceleration
    case 'm/s^2'
      %x = x;
    case 'cm/s^2'
      x = x.*100;
    case 'mm/s^2'
      x = x.*1000;
    case 'ft/s^2'
      x = x./0.3048;
    case 'in/s^2'
      x = x./0.0254;
    case 'G'
      x = x./9.81;

    % Area
    case 'km^2'
      x = x./1e6;
    case 'm^2'
      %x = x;
    case 'cm^2'
      x = x.*10^4;
    case 'mm^2'
      x = x.*10^6;
    case 'ym^2'
      x = x.*10^12;
    case 'sq-mile'
      x = x./1609.344^2;
    case 'sq-yd'
      x = x./0.83612736;
    case 'sq-ft'
      x = x./0.3048^2;
    case 'sq-in'
      x = x./0.0254^2;
    case 'acres'
      x = x./4046.8564224;
    case 'ha'
      x = x./1e4;
    case 'ar'
      x = x./100;

    % Volume
    case 'km^3'
      x = x./10^9;
    case 'm^3'
      %x = x;
    case 'cm^3'
      x = x.*10^6;
    case 'mm^3'
      x = x.*10^9;
    case 'ym^3'
      x = x.*10^18;
    case 'cu-mile'
      x = x./1609.344^3;
    case 'cu-ft'
      x = x./0.3048^3;
    case 'cu-yd'
      x = x./0.764554857984;
    case 'cu-in'
      x = x./0.0254^3;
    case 'l'
      x = x.*1000;
    case 'cl'
      x = x.*1e6;
    case 'ml'
      x = x.*10^9;
    case 'galUS'
      x = x./3.785412e-3;
    case 'galUK'
      x = x./4.54609e-3; 
    case 'pint'
      x = x./4.73176473e-4;
    case 'quart'
      x = x./9.46352946e-4;
    case 'fl-oz'
      x = x./2.95735295625e-5;

    % Density
    case 't/m^3'
      x = x./1000;
    case 'kg/m^3'
      %x = x;
    case 'g/cm^3'
      x = x./1000;
    case 'g/mm^3'
      x = x./1e6;
    case 'kg/l'
      x = x./1000;
    case 'lbs/ft^3'
      x = x./(0.4535924/0.3048^3);
    case 'lbs/in^3'
      x = x./(0.4535924/0.0254^3);
    case 'lbs/galUS'
      x = x./119.826427;
    case 'lbs/galUK'
      x = x./99.776373;

    % Mass Moment of Inertia
    case 'kg*m^2'
      %x = x;
    case 'kg*cm^2'
      x = x.*1e4;
    case 'kg*mm^2'
      x = x.*1e6;
    case 'g*m^2'
      x = x.*1000;
    case 'g*cm^2'
      x = x.*1e7;
    case 'g*mm^2'
      x = x.*1e9;
    case 'lb*ft^2'
      x = x./0.04214011;
    case 'lb*in^2'
      x = x./0.2926397e-3;

    % Areal Moment of Inertia
    case 'm^4'
      %x = x;
    case 'cm^4'
      x = x.*1e8;
    case 'mm^4'
      x = x.*1e12;
    case 'ft^4'
      x = x.*(0.3048^4);
    case 'in^4'
      x = x.*(0.0254^4);

    % Frequency / Angular Velocity
    case 'GHz'
      x = x./1e9;
    case 'MHz'
      x = x./1e6;
    case 'kHz'
      x = x./1e3;
    case 'Hz'
      %x = x;
    case '1/min'
      x = x.*60;
    case '1/h'
      x = x.*3600;
    case 'rad/s'
      x = x.*(2*pi);
    case 'deg/s'
      x = x.*360;
    case 'rpm'
      x = x.*60;

    % Angle
    case 'rad'
      %x = x;
    case 'deg'
      x = x.*(180/pi);
    case 'rev'
      x = x./(2*pi);

    % Stress Intensity Factor
    case 'MPa*m^1/2'
      %x = x;
    case 'MPa*mm^1/2'
      x = x.*sqrt(1000);
    case 'MPa*in^1/2'
      x = x./6.894757./sqrt(0.0254);

    % Fuel Consumption
    case 'l/100km'
      %x = x;
    case 'miles/galUS'
      x = 235.214596754951./x;

    % Viscosity
    case 'Ns/m^2'
      %x = x;
    case 'poise'
      x = x.*10;
    case 'centipoise'
      x = x.*1000;
    case 'lbfs/sq-ft'
      x = x.*0.02089;

    % Strain
    case 'm/m'
      %x = x;
    case 'mm/m'
      x = x.*1e3;
    case 'ym/m'
      x = x.*1e6;
    case 'nm/m'
      x = x.*1e9;
    case '%'
      x = x.*1e2;
    case 'o/oo'
      x = x.*1e3;
    case 'in/in'
      %x = x;
    case 'mill/in'
      x = x.*1e3';

    % Metric Prefixes
    case 'Yotta'
      x = x.*1e-24;
    case 'Zetta'
      x = x.*1e-21;
    case 'Exa'
      x = x.*1e-18;
    case 'Peta'
      x = x.*1e-15;
    case 'Tera'
      x = x.*1e-12;
    case 'Giga'
      x = x.*1e-9;
    case 'Mega'
      x = x.*1e-6;
    case 'Myria'
      x = x.*1e-5;
    case 'kilo'
      x = x.*1e-3;
    case 'hecto'
      x = x.*1e-2;
    case '1'
      %x = x;
    case 'deci'
      x = x.*1e1;
    case 'centi'
      x = x.*1e2;
    case 'milli'
      x = x.*1e3;
    case 'micro'
      x = x.*1e6;
    case 'nano'
      x = x.*1e9;
    case 'pico'
      x = x.*1e12;
    case 'femto'
      x = x.*1e15;
    case 'atto'
      x = x.*1e18;
    case 'zepto'
      x = x.*1e21;
    case 'yocto'
      x = x.*1e24;


otherwise
  write_error_to_log(['Error - Unsupported unit2: ', unit2], mfilename);

end
end

return