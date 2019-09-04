#ifndef UTILITIESTEST_H
#define UTILITIESTEST_H
#include "utilities.h"

#include "gtest/gtest.h"

#include "BasicTestClass.hpp"

namespace reseq{
	class utilitiesTest : public BasicTestClass{
	public: 
		static void Register();
	};
}

#endif // UTILITIESTEST_H
