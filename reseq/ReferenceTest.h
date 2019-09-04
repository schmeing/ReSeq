#ifndef REFERENCETEST_H
#define REFERENCETEST_H
#include "Reference.h"

#include <stdint.h>

#include "gtest/gtest.h"

#include "BasicTestClass.hpp"

namespace reseq{
	class ReferenceTest : public BasicTestClass{
	public:
		static void Register();

	protected:
		Reference ref_;

		void TestVariantClass();
		void TestInsertVariant();
		void TestVariationLoading();
		void TestVariationPositionLoading();
	};
}

#endif // REFERENCETEST_H
