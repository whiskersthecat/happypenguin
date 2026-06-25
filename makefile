happypenguin: happypenguin2.0.4.cpp happypenguin.h
	g++ happypenguin2.0.4.cpp -O3 -std=c++11 -Wall -Wextra -o happypenguin
debug: happypenguin2.0.4.cpp happypenguin.h
	g++ happypenguin2.0.4.cpp -g -O3 -std=c++11 -Wall -Wextra -o happypenguin
recombigator: recombigator.cc recombigator.h
	g++ recombigator.cc -O3 -std=c++11 -Wall -Wextra -o recombigator
checkSNPs: checkSNPs.cc
	g++ checkSNPs.cc -O3 -std=c++11 -o checkSNPs
