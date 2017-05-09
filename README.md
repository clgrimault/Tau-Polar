# Tau-Polar

git clone https://github.com/cherepan/Tau-Polar.git

and run ./todo.pl for an instruction


In order to include a user class to the combined library, compile your project in UserCodes directory and accordingly  modify tauola/example Makefile to be able to use your library :
USERLIBS = -L./UserCodes -lUserLib
