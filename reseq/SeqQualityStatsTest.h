#ifndef SEQQUALITYSTATSTEST_H
#define SEQQUALITYSTATSTEST_H
#include "SeqQualityStats.hpp"

#include "gtest/gtest.h"

#include "BasicTestClass.hpp"

namespace reseq{
	class SeqQualityStatsTest : public BasicTestClass{
	public: 
		static void Register();
	};
}

#endif // SEQQUALITYSTATSTEST_H
