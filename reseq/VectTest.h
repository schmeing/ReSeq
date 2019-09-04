#ifndef VECTTEST_H
#define VECTTEST_H
#include "Vect.hpp"

#include "gtest/gtest.h"

#include "BasicTestClass.hpp"

namespace reseq{
	class VectTest : public BasicTestClass{
	public: 
		static void Register();
	};
}

#endif // VECTTEST_H
