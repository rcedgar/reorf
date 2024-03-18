#include "myutils.h"
#include "alpha.h"

const char *KmerCodeToStr(uint KmerCode, uint k, string &Str)
	{
	Str.clear();
	for (uint i = 0; i < k; ++i)
		{
		uint Letter = KmerCode%4;
		Str += g_LetterToCharNucleo[Letter];
		KmerCode /= 4;
		}
	reverse(Str.begin(), Str.end());
	return Str.c_str();
	}

uint32 StrToKmerCode(const string &Seq)
	{
	uint k = SIZE(Seq);
	uint32 KmerCode = 0;
	for (uint i = 0; i < k; ++i)
		{
		char c = Seq[i];
		uint Letter = g_CharToLetterNucleo[c];
		if (Letter == INVALID_LETTER)
			return UINT32_MAX;
		KmerCode = (KmerCode << 2) | Letter;
		}
	return KmerCode;
	}

uint32 SeqToKmerCode(const byte *Seq, uint k)
	{
	uint32 KmerCode = 0;
	for (uint i = 0; i < k; ++i)
		{
		char c = Seq[i];
		uint Letter = g_CharToLetterNucleo[c];
		if (Letter == INVALID_LETTER)
			return UINT32_MAX;
		KmerCode = (KmerCode << 2) | Letter;
		}
	return KmerCode;
	}

void TestKmerCode(uint k)
	{
	const uint N = myipow(k, 4);
	for (uint KmerCode = 0; KmerCode < N; ++KmerCode)
		{
		string Str;
		KmerCodeToStr(KmerCode, k, Str);

		uint KmerCode2 = StrToKmerCode(Str);
		asserta(KmerCode2 == KmerCode);
		}
	ProgressLog("TestKmerCode ok\n");
	}
