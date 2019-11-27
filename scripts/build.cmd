::==========================================================================
:: Author: Oriol Cervelló (oriol.cn [at] protonmail.com) 
::==========================================================================
:: License: GNU GPLv3.0
:: Copyright (C) 2019  Oriol Cervelló
::
:: This program is free software: you can redistribute it and/or modify
:: it under the terms of the GNU General Public License as published by
:: the Free Software Foundation, either version 3 of the License, or
:: (at your option) any later version.
:: 
:: This program is distributed in the hope that it will be useful,
:: but WITHOUT ANY WARRANTY; without even the implied warranty of
:: MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
:: GNU General Public License for more details.
:: 
:: You should have received a copy of the GNU General Public License
:: along with this program.  If not, see <http://www.gnu.org/licenses/>.
::==========================================================================

IF %1 EQU 1 (
 nvcc ..\others\DeviceProperties.cu   -o ..\Others\DeviceProperties
 .\..\others\DeviceProperties
 exit /b
)

IF  %1 EQU 2 (
 nvcc ../src/TextParser.cu ../src/GlobalFunc.cu ../src/IOFunc.cu ../src/HostFunc.cu ../src/main.cu -lnpps -lcufft -o ../bin/gnssprocessgpu
 exit /b
)

IF  %1 EQU 3 (
 nvcc ../src/TextParser.cu ../src/GlobalFunc.cu ../src/IOFunc.cu ../src/HostFunc.cu ../src/mainFloat.cu -lnpps -lcufft -o ../bin/gnssprocessgpu
 exit /b
)

IF  %1 EQU 4 (
 nvcc ../src/TextParser.cu ../src/GlobalFunc.cu ../src/IOFunc.cu ../src/HostFunc.cu ../src/mainInt16.cu -lnpps -lcufft -o ../bin/gnssprocessgpu
 exit /b
)
