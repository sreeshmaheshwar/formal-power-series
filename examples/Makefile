CXX = g++
CXXFLAGS = -std=c++20 -Wall -Wextra -Wpedantic
INCLUDES = -I ../ac-library

single: $(file)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(file) -o $(basename $(file)).out

clean:
	rm -f *.out
	rm -f */*.out
