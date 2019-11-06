#ifndef SURROUNDINGTEST_H
#define SURROUNDINGTEST_H
#include "Surrounding.h"

#include "gtest/gtest.h"

#include "BasicTestClass.hpp"

namespace reseq{
	class SurroundingTest : public BasicTestClassWithReference{
	public: 
		static void Register();

	protected:
		virtual void TearDown();

		void TestBasics();
		void TestSettersAndUpdaters();
		void TestSettersAndUpdatersWithN();
		void TestModifiers();
		void TestModifiersExtremCases();
		void TestCombiningBias();
		void TestSeparatingBias();
	};
}

#endif // SURROUNDINGTEST_H
