#include "myutils.h"
#include "seqdb.h"
#include "translator.h"
#include "bipal.h"
#include "alpha.h"

uint GetCodonDiffCount(const char CodonA[3], const char CodonB[3])
	{
	uint n = 0;
	for (uint i = 0; i < 3; ++i)
		{
		if (CodonA[i] != CodonB[i])
			++n;
		}
	return n;
	}

void MakeAnnotRow(const string &Row1, const string &Row2, string &AnnotRow)
	{
	AnnotRow.clear();
	const uint ColCount = SIZE(Row1);
	asserta(SIZE(Row2) == ColCount);
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		if (Row1[Col] == Row2[Col])
			AnnotRow += '|';
		else
			AnnotRow += ' ';
		}
	}

void GetPosToCol(const string &Row, vector<uint> &PosToCol)
	{
	PosToCol.clear();
	const uint ColCount = SIZE(Row);
	uint Pos = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		if (Row[Col] != '-')
			PosToCol.push_back(Col);
		}
	}

void GetColToPos(const string &Row, vector<uint> &ColToPos)
	{
	ColToPos.clear();
	const uint ColCount = SIZE(Row);
	uint Pos = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		if (Row[Col] == '-')
			ColToPos.push_back(UINT_MAX);
		else
			ColToPos.push_back(Pos++);
		}
	}

void ValidateColToPos(const vector<uint> &ColToPos, uint L)
	{
	const uint ColCount = SIZE(ColToPos);

	uint PrevPos = UINT_MAX;
	vector<bool> PosFound(L, false);
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		uint Pos = ColToPos[Col];
		if (Pos == UINT_MAX)
			continue;
		asserta(Pos < L);
		if (PrevPos != UINT_MAX)
			asserta(Pos > PrevPos);
		PosFound[Pos] = true;
		PrevPos = Pos;
		}
	for (uint Pos = 0; Pos < L; ++Pos)
		asserta(PosFound[Pos]);
	}

void ValidateColPos(const vector<uint> &PosToCol, const vector<uint> &ColToPos)
	{
	const uint ColCount = SIZE(ColToPos);
	const uint L = SIZE(PosToCol);
	uint PrevCol = UINT_MAX;
	for (uint Pos = 0; Pos < L; ++Pos)
		{
		uint Col = PosToCol[Pos];
		asserta(Col < ColCount);
		asserta(Pos == ColToPos[Col]);
		if (Pos > 0)
			asserta(Col > PrevCol);
		PrevCol = Col;
		}

	uint PrevPos = UINT_MAX;
	vector<bool> PosFound(L, false);
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		uint Pos = ColToPos[Col];
		if (Pos == UINT_MAX)
			continue;
		asserta(Col == PosToCol[Pos]);
		if (PrevPos != UINT_MAX)
			asserta(Pos > PrevPos);
		PosFound[Pos] = true;
		PrevPos = Pos;
		}
	for (uint Pos = 0; Pos < L; ++Pos)
		asserta(PosFound[Pos]);
	}

void StripGaps(const string &Row, string &Seq)
	{
	Seq.clear();
	for (uint i = 0; i < SIZE(Row); ++i)
		{
		char c = Row[i];
		if (c != '-')
			Seq += c;
		}
	}

static void AssertNoStops(const string &AaStr)
	{
	for (uint i = 0; i < SIZE(AaStr); ++i)
		asserta(AaStr[i] != '*');
	}

void BiPal::InitFrames()
	{
	m_PlusFrame = 0;
	m_MinusFrame = 0;

	vector<int> EndlessFramesX;
	vector<int> EndlessFramesY;
	vector<bool> HasMetsX;
	vector<bool> HasMetsY;

	m_Translator.SetStr(m_NtPlusSeqX, true);
	m_Translator.GetEndlessFrames(EndlessFramesX, HasMetsX);
	m_Translator.GetEndlessFrames(EndlessFramesY, HasMetsY);

	for (uint i = 0; i < SIZE(EndlessFramesX); ++i)
		{
		int Frame = EndlessFramesX[i];
		if (Frame > 0)
			{
			asserta(m_PlusFrame == 0);
			m_PlusFrame = Frame;
			m_PlusFrameOffset = uint(Frame - 1);
			asserta(m_PlusFrameOffset <= 2);
			}
		else if (Frame < 0)
			{
			asserta(m_MinusFrame == 0);
			m_MinusFrame = Frame;
			m_MinusFrameOffset = uint(-Frame - 1);
			asserta(m_MinusFrameOffset <= 2);
			}
		else
			asserta(false);
		}

	for (uint i = 0; i < SIZE(EndlessFramesY); ++i)
		{
		int Frame = EndlessFramesY[i];
		if (Frame > 0)
			asserta(Frame == m_PlusFrame);
		else if (Frame < 0)
			asserta(Frame == m_MinusFrame);
		else
			asserta(false);
		}
	ProgressLog("Frames %+d, %+d (phase %d)\n",
	  m_PlusFrame, m_MinusFrame, GetPhase(m_PlusFrame, m_MinusFrame));
	}

void BiPal::InitOrfs()
	{
	m_NtPlusOrfX.clear();
	m_NtPlusOrfY.clear();

	m_NtMinusOrfX.clear();
	m_NtMinusOrfY.clear();

	int LX = int(m_LX);
	int LY = int(m_LY);
	for (uint i = 0; i < m_LX; ++i)
		{
		uint Pos = i + m_PlusFrameOffset;
		char c = m_NtPlusSeqX[Pos%m_LX];
		m_NtPlusOrfX += c;
		}

	for (int i = 0; i < LY; ++i)
		{
		uint Pos = i + m_PlusFrameOffset;
		char c = m_NtPlusSeqY[Pos%m_LY];
		m_NtPlusOrfY += c;
		}

	for (int i = 0; i < LX; ++i)
		{
		int iPos = LX - int(m_MinusFrameOffset) - i - 1;
		if (iPos < 0)
			iPos += LX;
		asserta(iPos >= 0 && iPos < LY);
		char c = m_NtPlusSeqX[iPos];
		char cc = g_CharToCompChar[c];
		m_NtMinusOrfX += cc;
		}

	for (int i = 0; i < LY; ++i)
		{
		int iPos = LY - int(m_MinusFrameOffset) - i - 1;
		if (iPos < 0)
			iPos += LY;
		asserta(iPos >= 0 && iPos < LY);
		char c = m_NtPlusSeqY[iPos];
		char cc = g_CharToCompChar[c];
		m_NtMinusOrfY += cc;
		}

	m_Translator.SetStr(m_NtPlusOrfX, false);
	m_Translator.SetFrame(1);
	m_Translator.GetAaStr(m_AaPlusOrfX);

	m_Translator.SetStr(m_NtPlusOrfY, false);
	m_Translator.SetFrame(1);
	m_Translator.GetAaStr(m_AaPlusOrfY);

	m_Translator.SetStr(m_NtMinusOrfX, false);
	m_Translator.SetFrame(1);
	m_Translator.GetAaStr(m_AaMinusOrfX);

	m_Translator.SetStr(m_NtMinusOrfY, false);
	m_Translator.SetFrame(1);
	m_Translator.GetAaStr(m_AaMinusOrfY);

	AssertNoStops(m_AaPlusOrfX);
	AssertNoStops(m_AaPlusOrfY);
	AssertNoStops(m_AaMinusOrfX);
	AssertNoStops(m_AaMinusOrfY);
	}

void BiPal::InitNtPlusPosToCol()
	{
	GetPosToCol(m_NtPlusRowX, m_NtPlusPosToColX);
	GetPosToCol(m_NtPlusRowY, m_NtPlusPosToColY);
	GetColToPos(m_NtPlusRowX, m_NtPlusColToPosX);
	GetColToPos(m_NtPlusRowY, m_NtPlusColToPosY);
	}

uint BiPal::PlusNtPosToMinusOrfAaPosX(uint PlusNtPosX, uint &BaseNr) const
	{
	asserta(PlusNtPosX < m_LX);
	uint MinusNtPosX = m_LX - PlusNtPosX - 1;
	uint MinusOrfNtPosX = (m_LX + MinusNtPosX - m_MinusFrameOffset)%m_LX;
	uint MinusOrfAaPosX = MinusOrfNtPosX/3;
	BaseNr = MinusOrfNtPosX%3;
	return MinusOrfAaPosX;
	}

uint BiPal::PlusNtPosToMinusOrfAaPosY(uint PlusNtPosY, uint &BaseNr) const
	{
	asserta(PlusNtPosY < m_LY);
	uint MinusNtPosY = m_LY - PlusNtPosY - 1;
	uint MinusOrfNtPosY = (m_LY + MinusNtPosY - m_MinusFrameOffset)%m_LY;
	uint MinusOrfAaPosY = MinusOrfNtPosY/3;
	BaseNr = MinusOrfNtPosY%3;
	return MinusOrfAaPosY;
	}

uint BiPal::PlusOrfAaPosToPlusNtPosX(uint PlusOrfAaPos, uint BaseNr) const
	{
	uint PlusOrfNtPosX = (PlusOrfAaPos*3 + BaseNr)%m_LX;
	uint PlusNtPosX = (PlusOrfNtPosX + m_LX - m_PlusFrameOffset)%m_LX;
	return PlusNtPosX;
	}

uint BiPal::PlusOrfAaPosToPlusNtPosY(uint PlusOrfAaPos, uint BaseNr) const
	{
	uint PlusOrfNtPosY = (PlusOrfAaPos*3 + BaseNr)%m_LY;
	uint PlusNtPosY = (PlusOrfNtPosY + m_LY - m_PlusFrameOffset)%m_LY;
	return PlusNtPosY;
	}

uint BiPal::MinusOrfAaPosToPlusNtPosX(uint MinusOrfAaPos, uint BaseNr) const
	{
	uint PlusOrfNtPosX = (2*m_LX - (3*MinusOrfAaPos + BaseNr) - 1)%m_LX;
	uint PlusNtPosX = (m_LX + PlusOrfNtPosX - m_MinusFrameOffset)%m_LX;
	return PlusNtPosX;
	}

uint BiPal::MinusOrfAaPosToPlusNtPosY(uint MinusOrfAaPos, uint BaseNr) const
	{
	uint PlusOrfNtPosY = (2*m_LY - (3*MinusOrfAaPos + BaseNr) - 1)%m_LY;
	uint PlusNtPosY = (m_LY + PlusOrfNtPosY - m_MinusFrameOffset)%m_LY;
	return PlusNtPosY;
	return 0;
	}

uint BiPal::PlusNtPosToPlusOrfAaPosX(uint PlusOrfNtPosX, uint &BaseNr) const
	{
	uint PlusNtPosX = (PlusOrfNtPosX + m_PlusFrameOffset)%m_LX;
	uint PlusOrfAaPosX = PlusNtPosX/3;
	BaseNr = PlusNtPosX%3;
	return PlusOrfAaPosX;
	}

uint BiPal::PlusNtPosToPlusOrfAaPosY(uint PlusOrfNtPosY, uint &BaseNr) const
	{
	uint PlusNtPosY = (PlusOrfNtPosY + m_PlusFrameOffset)%m_LY;
	uint PlusOrfAaPosY = PlusNtPosY/3;
	BaseNr = PlusNtPosY%3;
	return PlusOrfAaPosY;
	}

uint BiPal::CalcAaOrfLengthX() const
	{
	asserta(m_LX%3 == 0);
	return m_LX/3;
	}

uint BiPal::CalcAaOrfLengthY() const
	{
	asserta(m_LY%3 == 0);
	return m_LY/3;
	}

void BiPal::GetBarStr(const char *PlusSeq, uint Offset, string &BarStr) const
	{
	BarStr.clear();
	const uint L = ustrlen(PlusSeq);

	BarStr.clear();
	asserta(Offset <= 2);
	for (uint i = 0; i < 2 - Offset; ++i)
		BarStr += ' ';
	for (uint i = 0; i < Offset; ++i)
		BarStr += PlusSeq[i];

	BarStr += " |";
	BarStr += PlusSeq[Offset];
	BarStr += "| ";

	for (uint i = Offset + 1; i < L; ++i)
		BarStr += PlusSeq[i];

	while(SIZE(BarStr) < 9)
		BarStr +=  ' ';
	}

void BiPal::GetCodonPlusStr(const char *PlusSeq, uint Start, string &Str) const
	{
	Str.clear();

	char Codon[3];
	const uint L = ustrlen(PlusSeq);
	asserta(Start <= 2);
	for (uint i = 0; i < 2 - Start; ++i)
		Str += ' ';
	for (uint i = 0; i < Start; ++i)
		Str += tolower(PlusSeq[i]);
	Str += "[";
	for (uint i = 0; i < 3; ++i)
		{
		char c = PlusSeq[Start+i];
		Str += c;
		Codon[i] = c;
		}
	char aa = TranslateCodon3(Codon);
	Str += "=";
	Str += aa;
	Str += "]";
	for (uint i = Start+3; i < L; ++i)
		Str += tolower(PlusSeq[i]);
	while (SIZE(Str) < 12)
		Str += ' ';
	}

void BiPal::GetCodonMinusStr(const char *PlusSeq, uint Start, string &Str) const
	{
	Str.clear();

	char Codon[3];

	const uint L = ustrlen(PlusSeq);
	asserta(Start <= 2);
	for (uint i = 0; i < 2 - Start; ++i)
		Str += ' ';
	for (uint i = 0; i < Start; ++i)
		Str += tolower(PlusSeq[i]);
	Str += "[";
	for (uint i = 0; i < 3; ++i)
		{
		char c = PlusSeq[Start+i];
		Codon[2-i] = g_CharToCompChar[c];
		Str += tolower(c);
		}
	Str += "<";
	for (uint i = 0; i < 3; ++i)
		Str += Codon[i];
	char aa = TranslateCodon3(Codon);
	Str += "=";
	Str += aa;
	Str += "]";
	for (uint i = Start+3; i < L; ++i)
		Str += tolower(PlusSeq[i]);
	while (SIZE(Str) < 15)
		Str += ' ';
	}

void BiPal::WriteCodonInfo(FILE *f, const char *PlusSeq, uint Offset,
  uint StartPlus, uint StartMinus) const
	{
	if (f == 0)
		return;

	string BarStr;
	GetBarStr(PlusSeq, Offset, BarStr);
	fprintf(f, "%s", BarStr.c_str());

	string CodonPlusStr;
	GetCodonPlusStr(PlusSeq, StartPlus, CodonPlusStr);
	fprintf(f, "  %s", CodonPlusStr.c_str());

	GetCodonMinusStr(PlusSeq, StartMinus, CodonPlusStr);
	fprintf(f, "  %s", CodonPlusStr.c_str());
	}

void BiPal::WriteCodons(FILE *f) const
	{
	if (f == 0)
		return;

	uint NtPlusPosX = 0;
	uint NtPlusPosY = 0;
	const uint ColCountNt = SIZE(m_NtPlusRowX);
	asserta(SIZE(m_NtPlusRowY) == ColCountNt);
	fprintf(f,
"\n"
"nt aln                        X                                           Y                          plus codons        minus codons\n"
"------       -------------------------------------       -------------------------------------      --------------     --------------\n"
"X   Y        Base     +codon aa+        -codon aa-       Base     +codon aa+        -codon aa-      xy  XXX YYY xy     yy  XXX YYY xy\n"
"|   |          |         ||| |             ||| |           |         ||| |             ||| |        ||  ||| ||| ||     ||  ||| ||| ||\n"
"\n"
);

	for (uint Col = 0; Col < ColCountNt; ++Col)
		{
		char x = m_NtPlusRowX[Col];
		char y = m_NtPlusRowY[Col];
		char AnnotChar = ' ';
		if (x != '-' && y != '-')
			{
			if (x == y)
				AnnotChar = '=';
			else
				AnnotChar = '.';
			}

		fprintf(f, "%c %c %c  ", x, AnnotChar, y);

		if (x != '-' && y != '-')
			{
			char PlusSeqX[6];
			char PlusSeqY[6];

			uint OffsetX = UINT_MAX;
			uint OffsetY = UINT_MAX;

			uint StartPlusX = UINT_MAX;
			uint StartPlusY = UINT_MAX;

			uint StartMinusX = UINT_MAX;
			uint StartMinusY = UINT_MAX;

			char AaPlusX, AaMinusX;
			char AaPlusY, AaMinusY;

			char CodonPlusX[3];
			char CodonPlusY[3];
			char CodonMinusX[3];
			char CodonMinusY[3];

			uint BaseNrPlusX = GetCodonPlusX(NtPlusPosX, CodonPlusX);
			uint BaseNrPlusY = GetCodonPlusY(NtPlusPosY, CodonPlusY);

			uint BaseNrMinusX = GetCodonMinusX(NtPlusPosX, CodonMinusX);
			uint BaseNrMinusY = GetCodonMinusY(NtPlusPosY, CodonMinusY);

			GetCodonInfoPlusMinusX(NtPlusPosX, PlusSeqX, AaPlusX, AaMinusX,
			  OffsetX, StartPlusX, StartMinusX);

			fprintf(f, "    ");
			WriteCodonInfo(g_fLog, PlusSeqX, OffsetX, StartPlusX, StartMinusX);

			fprintf(f, "    ");
			GetCodonInfoPlusMinusY(NtPlusPosY, PlusSeqY, AaPlusY, AaMinusY,
			  OffsetY, StartPlusY, StartMinusY);
			WriteCodonInfo(g_fLog, PlusSeqY, OffsetY, StartPlusY, StartMinusY);

			fprintf(f, "    ");
			fprintf(f, " %u", BaseNrPlusX);
			fprintf(f, "%u", BaseNrPlusY);
			fprintf(f, "  %3.3s %3.3s ", CodonPlusX, CodonPlusY);
			fprintf(f, "%c", BaseNrPlusX == 0 ? AaPlusX : ' ');
			fprintf(f, "%c", BaseNrPlusY == 0 ? AaPlusY : ' ');

			fprintf(f, "    ");
			fprintf(f, " %u", BaseNrMinusX);
			fprintf(f, "%u", BaseNrMinusY);
			fprintf(f, "  %3.3s %3.3s ", CodonMinusX, CodonMinusY);
			fprintf(f, "%c", BaseNrMinusX == 0 ? AaMinusX : ' ');
			fprintf(f, "%c", BaseNrMinusY == 0 ? AaMinusY : ' ');
			}

		if (x != '-')
			++NtPlusPosX;
		if (y != '-')
			++NtPlusPosY;

		fprintf(f, "\n");
		}
	}

void BiPal::LogMaps() const
	{
	const uint AaOrfLengthX = CalcAaOrfLengthX();
	const uint AaOrfLengthY = CalcAaOrfLengthY();

	uint NtPlusPosX = 0;
	uint NtPlusPosY = 0;
	const uint ColCount = SIZE(m_NtPlusRowX);
	asserta(SIZE(m_NtPlusRowY) == ColCount);
	uint BaseNr = UINT_MAX;

	Log(" Col   NtX   NtY  x y  +aaX/b  -aaX/b  +aaY/b  -aaY/b\n");
	//   1234  1234  1234  1 1  1234/1  1234/1  1234/1  1234/1
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char x = m_NtPlusRowX[Col];
		char y = m_NtPlusRowY[Col];

		uint AaPlusOrfPosX = UINT_MAX;
		uint AaMinusOrfPosX = UINT_MAX;
		Log("%4u  %4u  %4u  %c %c", Col, NtPlusPosX, NtPlusPosY, x, y);
		if (x != '-')
			{
			AaPlusOrfPosX = PlusNtPosToPlusOrfAaPosX(NtPlusPosX, BaseNr);
			Log("  %4u/%u", AaPlusOrfPosX, BaseNr);

			AaMinusOrfPosX = PlusNtPosToMinusOrfAaPosX(NtPlusPosX, BaseNr);
			Log("  %4u/%u", AaMinusOrfPosX, BaseNr);
			}
		else
			{
			Log("        ");
			Log("        ");
			}

		if (y != '-')
			{
			uint AaPlusOrfPosY = PlusNtPosToPlusOrfAaPosY(NtPlusPosY, BaseNr);
			Log("  %4u/%u", AaPlusOrfPosY, BaseNr);

			uint AaMinusOrfPosY = PlusNtPosToMinusOrfAaPosY(NtPlusPosY, BaseNr);
			Log("  %4u/%u", AaMinusOrfPosY, BaseNr);
			}
		else
			{
			Log("        ");
			Log("        ");
			}

		if (AaPlusOrfPosX != UINT_MAX)
			{
			uint PosY = m_AaPlusOrfPosXToPosY[AaPlusOrfPosX];
			if (PosY != UINT_MAX)
				Log(" +aaORF(X%u:Y%u)", AaPlusOrfPosX, PosY);
			}
		if (AaMinusOrfPosX != UINT_MAX)
			{
			uint PosY = m_AaMinusOrfPosXToPosY[AaMinusOrfPosX];
			if (PosY != UINT_MAX)
				Log(" -aaORF(X%u:Y%u)", AaMinusOrfPosX, PosY);
			}

		if (x != '-')
			++NtPlusPosX;
		if (y != '-')
			++NtPlusPosY;

		Log("\n");
		}
	asserta(NtPlusPosX == m_LX);
	asserta(NtPlusPosY == m_LY);
	}

void BiPal::InitAaOrfPosToPosPlus()
	{
	const uint AaOrfLengthX = CalcAaOrfLengthX();
	const uint AaOrfLengthY = CalcAaOrfLengthY();

	vector<vector<uint> > AaPlusOrfPosXToPosY;
	AaPlusOrfPosXToPosY.resize(AaOrfLengthX);
	for (uint i = 0; i < AaOrfLengthX; ++i)
		AaPlusOrfPosXToPosY[i].resize(3, UINT_MAX);

	uint NtPlusPosX = 0;
	uint NtPlusPosY = 0;
	const uint ColCount = SIZE(m_NtPlusRowX);
	asserta(SIZE(m_NtPlusRowY) == ColCount);
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char x = m_NtPlusRowX[Col];
		char y = m_NtPlusRowY[Col];

		if (x != '-' && y != '-')
			{
			uint BaseNrX = UINT_MAX;
			uint BaseNrY = UINT_MAX;
			uint AaPlusOrfPosX = PlusNtPosToPlusOrfAaPosX(NtPlusPosX, BaseNrX);
			uint AaPlusOrfPosY = PlusNtPosToPlusOrfAaPosY(NtPlusPosY, BaseNrY);

			asserta(AaPlusOrfPosX < AaOrfLengthX);
			asserta(AaPlusOrfPosY < AaOrfLengthY);
			asserta(BaseNrX < 3);
			asserta(BaseNrY < 3);

			asserta(AaPlusOrfPosXToPosY[AaPlusOrfPosX][BaseNrX] == UINT_MAX);
			AaPlusOrfPosXToPosY[AaPlusOrfPosX][BaseNrX] = AaPlusOrfPosY;
			}

		if (x != '-')
			++NtPlusPosX;
		if (y != '-')
			++NtPlusPosY;
		}
	asserta(NtPlusPosX == m_LX);
	asserta(NtPlusPosY == m_LY);

	m_AaPlusOrfPosXToPosY.clear();
	for (uint AaPlusOrfPosX = 0; AaPlusOrfPosX < AaOrfLengthX; ++AaPlusOrfPosX)
		{
		uint PosY0 = AaPlusOrfPosXToPosY[AaPlusOrfPosX][0];
		uint PosY1 = AaPlusOrfPosXToPosY[AaPlusOrfPosX][1];
		uint PosY2 = AaPlusOrfPosXToPosY[AaPlusOrfPosX][2];
		if (PosY0 == PosY1 && PosY0 == PosY2 && PosY0 != UINT_MAX)
			m_AaPlusOrfPosXToPosY.push_back(PosY0);
		else
			m_AaPlusOrfPosXToPosY.push_back(UINT_MAX);
		}
	}

void BiPal::InitAaOrfPosToPosMinus()
	{
	const uint AaOrfLengthX = CalcAaOrfLengthX();
	const uint AaOrfLengthY = CalcAaOrfLengthY();

	vector<vector<uint> > AaMinusOrfPosXToPosY;
	AaMinusOrfPosXToPosY.resize(AaOrfLengthX);
	for (uint i = 0; i < AaOrfLengthX; ++i)
		AaMinusOrfPosXToPosY[i].resize(3, UINT_MAX);

	uint NtPlusPosX = 0;
	uint NtPlusPosY = 0;
	const uint ColCount = SIZE(m_NtPlusRowX);
	asserta(SIZE(m_NtPlusRowY) == ColCount);
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char x = m_NtPlusRowX[Col];
		char y = m_NtPlusRowY[Col];

		if (x != '-' && y != '-')
			{
			uint BaseNrX = UINT_MAX;
			uint BaseNrY = UINT_MAX;
			uint AaMinusOrfPosX = PlusNtPosToMinusOrfAaPosX(NtPlusPosX, BaseNrX);
			uint AaMinusOrfPosY = PlusNtPosToMinusOrfAaPosY(NtPlusPosY, BaseNrY);

			asserta(AaMinusOrfPosX < AaOrfLengthX);
			asserta(AaMinusOrfPosY < AaOrfLengthY);
			asserta(BaseNrX < 3);
			asserta(BaseNrY < 3);

			asserta(AaMinusOrfPosXToPosY[AaMinusOrfPosX][BaseNrX] == UINT_MAX);
			AaMinusOrfPosXToPosY[AaMinusOrfPosX][BaseNrX] = AaMinusOrfPosY;
			}

		if (x != '-')
			++NtPlusPosX;
		if (y != '-')
			++NtPlusPosY;
		}
	asserta(NtPlusPosX == m_LX);
	asserta(NtPlusPosY == m_LY);

	m_AaMinusOrfPosXToPosY.clear();
	for (uint AaMinusOrfPosX = 0; AaMinusOrfPosX < AaOrfLengthX; ++AaMinusOrfPosX)
		{
		uint PosY0 = AaMinusOrfPosXToPosY[AaMinusOrfPosX][0];
		uint PosY1 = AaMinusOrfPosXToPosY[AaMinusOrfPosX][1];
		uint PosY2 = AaMinusOrfPosXToPosY[AaMinusOrfPosX][2];
		if (PosY0 == PosY1 && PosY0 == PosY2 && PosY0 != UINT_MAX)
			m_AaMinusOrfPosXToPosY.push_back(PosY0);
		else
			m_AaMinusOrfPosXToPosY.push_back(UINT_MAX);
		}
	}

void BiPal::InitAaRowsPlus()
	{
	m_AaPlusRowX.clear();
	m_AaPlusRowY.clear();

	const uint AaOrfLengthX = CalcAaOrfLengthX();
	const uint AaOrfLengthY = CalcAaOrfLengthY();

	asserta(SIZE(m_AaPlusOrfPosXToPosY) == AaOrfLengthX);

	int iLastPosX = -1;
	int iLastPosY = -1;
	for (uint PosX = 0; PosX < AaOrfLengthX; ++PosX)
		{
		int iPosX = int(PosX);
		uint PosY = m_AaPlusOrfPosXToPosY[PosX];
		if (PosY == UINT_MAX)
			{
			asserta(PosX < SIZE(m_AaPlusOrfX));
			char x = m_AaPlusOrfX[PosX];
			m_AaPlusRowX += x;
			m_AaPlusRowY += '-';
			iLastPosX = int(PosX);
			continue;
			}

		int iPosY = int(PosY);
		asserta(iPosY > iLastPosY);
		for (int i = iLastPosY + 1; i < iPosY; ++i)
			{
			asserta(uint(i) < SIZE(m_AaPlusOrfY));
			char y = m_AaPlusOrfY[i];
			m_AaPlusRowX += '-';
			m_AaPlusRowY += y;
			}

		asserta(PosX < SIZE(m_AaPlusOrfX));
		char x = m_AaPlusOrfX[PosX];
		asserta(PosY < SIZE(m_AaPlusOrfY));
		char y = m_AaPlusOrfY[PosY];
		m_AaPlusRowX += x;
		m_AaPlusRowY += y;

		iLastPosX = iPosX;
		iLastPosY = iPosY;
		}

	GetColToPos(m_AaPlusRowX, m_AaPlusColToPosX);
	GetColToPos(m_AaPlusRowY, m_AaPlusColToPosY);

	GetPosToCol(m_AaPlusRowX, m_AaPlusPosToColX);
	GetPosToCol(m_AaPlusRowY, m_AaPlusPosToColY);

	ValidateColPos(m_AaPlusPosToColX, m_AaPlusColToPosX);
	ValidateColPos(m_AaPlusPosToColY, m_AaPlusColToPosY);

	Log("%s\n", m_AaPlusRowX.c_str());
	Log("%s\n", m_AaPlusRowY.c_str());
	}

void BiPal::InitAaRowsMinus()
	{
	m_AaMinusRowX.clear();
	m_AaMinusRowY.clear();

	const uint AaOrfLengthX = CalcAaOrfLengthX();
	const uint AaOrfLengthY = CalcAaOrfLengthY();

	asserta(SIZE(m_AaMinusOrfPosXToPosY) == AaOrfLengthX);

	int iLastPosX = -1;
	int iLastPosY = -1;
	for (uint PosX = 0; PosX < AaOrfLengthX; ++PosX)
		{
		int iPosX = int(PosX);
		uint PosY = m_AaMinusOrfPosXToPosY[PosX];
		if (PosY == UINT_MAX)
			{
			asserta(PosX < SIZE(m_AaMinusOrfX));
			char x = m_AaMinusOrfX[PosX];
			m_AaMinusRowX += x;
			m_AaMinusRowY += '-';
			iLastPosX = int(PosX);
			continue;
			}

		int iPosY = int(PosY);
		asserta(iPosY > iLastPosY);
		for (int i = iLastPosY + 1; i < iPosY; ++i)
			{
			asserta(uint(i) < SIZE(m_AaMinusOrfY));
			char y = m_AaMinusOrfY[i];
			m_AaMinusRowX += '-';
			m_AaMinusRowY += y;
			}

		asserta(PosX < SIZE(m_AaMinusOrfX));
		char x = m_AaMinusOrfX[PosX];
		asserta(PosY < SIZE(m_AaMinusOrfY));
		char y = m_AaMinusOrfY[PosY];
		m_AaMinusRowX += x;
		m_AaMinusRowY += y;

		iLastPosX = iPosX;
		iLastPosY = iPosY;
		}

	GetColToPos(m_AaMinusRowX, m_AaMinusColToPosX);
	GetColToPos(m_AaMinusRowY, m_AaMinusColToPosY);

	GetPosToCol(m_AaMinusRowX, m_AaMinusPosToColX);
	GetPosToCol(m_AaMinusRowY, m_AaMinusPosToColY);

	ValidateColPos(m_AaMinusPosToColX, m_AaMinusColToPosX);
	ValidateColPos(m_AaMinusPosToColY, m_AaMinusColToPosY);

	Log("%s\n", m_AaMinusRowX.c_str());
	Log("%s\n", m_AaMinusRowY.c_str());
	}

void BiPal::ValidateMaps() const
	{
	uint AaOrfLengthX = CalcAaOrfLengthX();
	uint AaOrfLengthY = CalcAaOrfLengthY();
	
	for (uint AaOrfPosX = 0; AaOrfPosX < AaOrfLengthX; ++AaOrfPosX)
		{
		for (uint BaseNr = 0; BaseNr < 3; ++BaseNr)
			{
			uint NtPosX = PlusOrfAaPosToPlusNtPosX(AaOrfPosX, BaseNr);

			uint BaseNr2 = UINT_MAX;
			uint AaOrfPosX2 = PlusNtPosToPlusOrfAaPosX(NtPosX, BaseNr2);
			asserta(AaOrfPosX2 == AaOrfPosX);
			asserta(BaseNr2 == BaseNr);
			}
		}

	for (uint AaOrfPosY = 0; AaOrfPosY < AaOrfLengthY; ++AaOrfPosY)
		{
		for (uint BaseNr = 0; BaseNr < 3; ++BaseNr)
			{
			uint NtPosY = PlusOrfAaPosToPlusNtPosY(AaOrfPosY, BaseNr);

			uint BaseNr2 = UINT_MAX;
			uint AaOrfPosY2 = PlusNtPosToPlusOrfAaPosY(NtPosY, BaseNr2);
			asserta(AaOrfPosY2 == AaOrfPosY);
			asserta(BaseNr2 == BaseNr);
			}
		}

	for (uint AaOrfPosX = 0; AaOrfPosX < AaOrfLengthX; ++AaOrfPosX)
		{
		for (uint BaseNr = 0; BaseNr < 3; ++BaseNr)
			{
			uint NtPosX = MinusOrfAaPosToPlusNtPosX(AaOrfPosX, BaseNr);

			uint BaseNr2 = UINT_MAX;
			uint AaOrfPosX2 = PlusNtPosToMinusOrfAaPosX(NtPosX, BaseNr2);
			asserta(AaOrfPosX2 == AaOrfPosX);
			asserta(BaseNr2 == BaseNr);
			}
		}

	for (uint AaOrfPosY = 0; AaOrfPosY < AaOrfLengthY; ++AaOrfPosY)
		{
		for (uint BaseNr = 0; BaseNr < 3; ++BaseNr)
			{
			uint NtPosY = MinusOrfAaPosToPlusNtPosY(AaOrfPosY, BaseNr);

			uint BaseNr2 = UINT_MAX;
			uint AaOrfPosY2 = PlusNtPosToMinusOrfAaPosY(NtPosY, BaseNr2);
			asserta(AaOrfPosY2 == AaOrfPosY);
			asserta(BaseNr2 == BaseNr);
			}
		}
	}

void BiPal::GetCodonInfoPlusMinusX(
  uint PlusNtPosX, char *PlusSeq,
  char &AaPlus, char &AaMinus,
  uint &Offset, uint &StartPlus, uint &StartMinus) const
	{
	uint BaseNrPlus = UINT_MAX;
	uint BaseNrMinus = UINT_MAX;

	uint AaOrfPosPlus = PlusNtPosToPlusOrfAaPosX(PlusNtPosX, BaseNrPlus);
	uint AaOrfPosMinus = PlusNtPosToMinusOrfAaPosX(PlusNtPosX, BaseNrMinus);
	
	asserta(AaOrfPosPlus < SIZE(m_AaPlusOrfX));
	AaPlus = m_AaPlusOrfX[AaOrfPosPlus];

	asserta(AaOrfPosMinus < SIZE(m_AaMinusOrfX));
	AaMinus = m_AaMinusOrfX[AaOrfPosMinus];

	uint PlusNtPosX_PlusBase0 = PlusOrfAaPosToPlusNtPosX(AaOrfPosPlus, 0);
	uint PlusNtPosX_PlusBase1 = PlusOrfAaPosToPlusNtPosX(AaOrfPosPlus, 1);
	uint PlusNtPosX_PlusBase2 = PlusOrfAaPosToPlusNtPosX(AaOrfPosPlus, 2);

	uint PlusNtPosX_MinusBase0 = MinusOrfAaPosToPlusNtPosX(AaOrfPosMinus, 0);
	uint PlusNtPosX_MinusBase1 = MinusOrfAaPosToPlusNtPosX(AaOrfPosMinus, 1);
	uint PlusNtPosX_MinusBase2 = MinusOrfAaPosToPlusNtPosX(AaOrfPosMinus, 2);

	uint Phase = GetPhase(m_PlusFrame, m_MinusFrame);
	asserta(Phase == 3);
	switch (BaseNrPlus)
		{
	case 0:
		asserta(PlusNtPosX_PlusBase0 == PlusNtPosX_MinusBase1);
		asserta(PlusNtPosX_PlusBase1 == PlusNtPosX_MinusBase0);
		PlusSeq[0] = m_NtPlusSeqX[PlusNtPosX_MinusBase2];
		PlusSeq[1] = m_NtPlusSeqX[PlusNtPosX_PlusBase0];
		PlusSeq[2] = m_NtPlusSeqX[PlusNtPosX_PlusBase1];
		PlusSeq[3] = m_NtPlusSeqX[PlusNtPosX_PlusBase2];
		PlusSeq[4] = 0;

		Offset = 1 + BaseNrPlus;
		StartPlus = 1;
		StartMinus = 0;
		break;

	case 1:
		asserta(PlusNtPosX_PlusBase1 == PlusNtPosX_MinusBase0);
		asserta(PlusNtPosX_PlusBase0 == PlusNtPosX_MinusBase1);
		PlusSeq[0] = m_NtPlusSeqX[PlusNtPosX_MinusBase2];
		PlusSeq[1] = m_NtPlusSeqX[PlusNtPosX_PlusBase0];
		PlusSeq[2] = m_NtPlusSeqX[PlusNtPosX_PlusBase1];
		PlusSeq[3] = m_NtPlusSeqX[PlusNtPosX_PlusBase2];
		PlusSeq[4] = 0;

		Offset = 1 + BaseNrPlus;
		StartPlus = 1;
		StartMinus = 0;
		break;

	case 2:
		asserta(PlusNtPosX_PlusBase2 == PlusNtPosX_MinusBase2);
		PlusSeq[0] = m_NtPlusSeqX[PlusNtPosX_PlusBase0];
		PlusSeq[1] = m_NtPlusSeqX[PlusNtPosX_PlusBase1];
		PlusSeq[2] = m_NtPlusSeqX[PlusNtPosX_PlusBase2];
		PlusSeq[3] = m_NtPlusSeqX[PlusNtPosX_MinusBase1];
		PlusSeq[4] = m_NtPlusSeqX[PlusNtPosX_MinusBase0];
		PlusSeq[5] = 0;

		Offset = BaseNrPlus;
		StartPlus = 0;
		StartMinus = 2;
		break;

	default:
		asserta(false);
		}

	char PlusCodon[3];
	for (uint i = 0; i < 3; ++i)
		{
		char c = PlusSeq[StartPlus+i];
		asserta(isalpha(c));
		PlusCodon[i] = c;
		}
	char AaPlus2 = TranslateCodon3(PlusCodon);

	char MinusCodon[3];
	for (uint i = 0; i < 3; ++i)
		{
		char c = PlusSeq[StartMinus+2-i];
		asserta(isalpha(c));
		MinusCodon[i] = g_CharToCompChar[c];
		}
	char AaMinus2 = TranslateCodon3(MinusCodon);

	asserta(AaPlus2 == AaPlus);
	asserta(AaMinus2 == AaMinus);
	}

void BiPal::GetCodonInfoPlusMinusY(
  uint PlusNtPosY, char *PlusSeq,
  char &AaPlus, char &AaMinus,
  uint &Offset, uint &StartPlus, uint &StartMinus) const
	{
	uint BaseNrPlus = UINT_MAX;
	uint BaseNrMinus = UINT_MAX;

	uint AaOrfPosPlus = PlusNtPosToPlusOrfAaPosY(PlusNtPosY, BaseNrPlus);
	uint AaOrfPosMinus = PlusNtPosToMinusOrfAaPosY(PlusNtPosY, BaseNrMinus);
	
	asserta(AaOrfPosPlus < SIZE(m_AaPlusOrfY));
	AaPlus = m_AaPlusOrfY[AaOrfPosPlus];

	asserta(AaOrfPosMinus < SIZE(m_AaMinusOrfY));
	AaMinus = m_AaMinusOrfY[AaOrfPosMinus];

	uint PlusNtPosY_PlusBase0 = PlusOrfAaPosToPlusNtPosY(AaOrfPosPlus, 0);
	uint PlusNtPosY_PlusBase1 = PlusOrfAaPosToPlusNtPosY(AaOrfPosPlus, 1);
	uint PlusNtPosY_PlusBase2 = PlusOrfAaPosToPlusNtPosY(AaOrfPosPlus, 2);

	uint PlusNtPosY_MinusBase0 = MinusOrfAaPosToPlusNtPosY(AaOrfPosMinus, 0);
	uint PlusNtPosY_MinusBase1 = MinusOrfAaPosToPlusNtPosY(AaOrfPosMinus, 1);
	uint PlusNtPosY_MinusBase2 = MinusOrfAaPosToPlusNtPosY(AaOrfPosMinus, 2);

	uint Phase = GetPhase(m_PlusFrame, m_MinusFrame);
	asserta(Phase == 3);
	switch (BaseNrPlus)
		{
	case 0:
		asserta(PlusNtPosY_PlusBase0 == PlusNtPosY_MinusBase1);
		asserta(PlusNtPosY_PlusBase1 == PlusNtPosY_MinusBase0);
		PlusSeq[0] = m_NtPlusSeqY[PlusNtPosY_MinusBase2];
		PlusSeq[1] = m_NtPlusSeqY[PlusNtPosY_PlusBase0];
		PlusSeq[2] = m_NtPlusSeqY[PlusNtPosY_PlusBase1];
		PlusSeq[3] = m_NtPlusSeqY[PlusNtPosY_PlusBase2];
		PlusSeq[4] = 0;

		Offset = 1 + BaseNrPlus;
		StartPlus = 1;
		StartMinus = 0;
		break;

	case 1:
		asserta(PlusNtPosY_PlusBase1 == PlusNtPosY_MinusBase0);
		asserta(PlusNtPosY_PlusBase0 == PlusNtPosY_MinusBase1);
		PlusSeq[0] = m_NtPlusSeqY[PlusNtPosY_MinusBase2];
		PlusSeq[1] = m_NtPlusSeqY[PlusNtPosY_PlusBase0];
		PlusSeq[2] = m_NtPlusSeqY[PlusNtPosY_PlusBase1];
		PlusSeq[3] = m_NtPlusSeqY[PlusNtPosY_PlusBase2];
		PlusSeq[4] = 0;

		Offset = 1 + BaseNrPlus;
		StartPlus = 1;
		StartMinus = 0;
		break;

	case 2:
		asserta(PlusNtPosY_PlusBase2 == PlusNtPosY_MinusBase2);
		PlusSeq[0] = m_NtPlusSeqY[PlusNtPosY_PlusBase0];
		PlusSeq[1] = m_NtPlusSeqY[PlusNtPosY_PlusBase1];
		PlusSeq[2] = m_NtPlusSeqY[PlusNtPosY_PlusBase2];
		PlusSeq[3] = m_NtPlusSeqY[PlusNtPosY_MinusBase1];
		PlusSeq[4] = m_NtPlusSeqY[PlusNtPosY_MinusBase0];
		PlusSeq[5] = 0;

		Offset = BaseNrPlus;
		StartPlus = 0;
		StartMinus = 2;
		break;

	default:
		asserta(false);
		}

	char PlusCodon[3];
	for (uint i = 0; i < 3; ++i)
		{
		char c = PlusSeq[StartPlus+i];
		asserta(isalpha(c));
		PlusCodon[i] = c;
		}
	char AaPlus2 = TranslateCodon3(PlusCodon);

	char MinusCodon[3];
	for (uint i = 0; i < 3; ++i)
		{
		char c = PlusSeq[StartMinus+2-i];
		asserta(isalpha(c));
		MinusCodon[i] = g_CharToCompChar[c];
		}
	char AaMinus2 = TranslateCodon3(MinusCodon);

	asserta(AaPlus2 == AaPlus);
	asserta(AaMinus2 == AaMinus);
	}

uint BiPal::GetCodonPlusX(uint PlusNtPos, char Codon[3]) const
	{
	asserta(PlusNtPos < m_LX);
	uint BaseNr = UINT_MAX;
	uint AaOrfPos = PlusNtPosToPlusOrfAaPosX(PlusNtPos, BaseNr);
	asserta(AaOrfPos < m_LX/3);

	uint CodonStartPos = (BaseNr > PlusNtPos ? (m_LX + PlusNtPos) - BaseNr : PlusNtPos - BaseNr);
	Codon[0] = m_NtPlusSeqX[CodonStartPos];
	Codon[1] = m_NtPlusSeqX[(CodonStartPos+1)%m_LX];
	Codon[2] = m_NtPlusSeqX[(CodonStartPos+2)%m_LX];

	return BaseNr;
	}

uint BiPal::GetCodonPlusY(uint PlusNtPos, char Codon[3]) const
	{
	asserta(PlusNtPos < m_LY);
	uint BaseNr = UINT_MAX;
	uint AaOrfPos = PlusNtPosToPlusOrfAaPosY(PlusNtPos, BaseNr);
	asserta(AaOrfPos < m_LY/3);

	uint CodonStartPos = (BaseNr > PlusNtPos ? (m_LY + PlusNtPos) - BaseNr : PlusNtPos - BaseNr);
	Codon[0] = m_NtPlusSeqY[CodonStartPos];
	Codon[1] = m_NtPlusSeqY[(CodonStartPos+1)%m_LY];
	Codon[2] = m_NtPlusSeqY[(CodonStartPos+2)%m_LY];

	return BaseNr;
	}

uint BiPal::GetCodonMinusX(uint PlusNtPos, char Codon[3]) const
	{
	asserta(PlusNtPos < m_LX);
	uint BaseNr = UINT_MAX;
	uint AaOrfPos = PlusNtPosToMinusOrfAaPosX(PlusNtPos, BaseNr);
	asserta(AaOrfPos < m_LX/3);
	asserta(AaOrfPos < SIZE(m_AaPlusOrfX));
	char aa = m_AaMinusOrfX[AaOrfPos];

	uint NtPos0 = MinusOrfAaPosToPlusNtPosX(AaOrfPos, 0);
	uint NtPos1 = MinusOrfAaPosToPlusNtPosX(AaOrfPos, 1);
	uint NtPos2 = MinusOrfAaPosToPlusNtPosX(AaOrfPos, 2);

	char b0 = m_NtPlusSeqX[NtPos0];
	char b1 = m_NtPlusSeqX[NtPos1];
	char b2 = m_NtPlusSeqX[NtPos2];

	Codon[0] = g_CharToCompChar[b0];
	Codon[1] = g_CharToCompChar[b1];
	Codon[2] = g_CharToCompChar[b2];

	char aa2 = TranslateCodon3(Codon);
	asserta(aa2 == aa);

	return BaseNr;
	}

uint BiPal::GetCodonMinusY(uint PlusNtPos, char Codon[3]) const
	{
	asserta(PlusNtPos < m_LY);
	uint BaseNr = UINT_MAX;
	uint AaOrfPos = PlusNtPosToMinusOrfAaPosY(PlusNtPos, BaseNr);
	asserta(AaOrfPos < m_LY/3);
	asserta(AaOrfPos < SIZE(m_AaPlusOrfY));
	char aa = m_AaMinusOrfY[AaOrfPos];

	uint NtPos0 = MinusOrfAaPosToPlusNtPosY(AaOrfPos, 0);
	uint NtPos1 = MinusOrfAaPosToPlusNtPosY(AaOrfPos, 1);
	uint NtPos2 = MinusOrfAaPosToPlusNtPosY(AaOrfPos, 2);

	char b0 = m_NtPlusSeqY[NtPos0];
	char b1 = m_NtPlusSeqY[NtPos1];
	char b2 = m_NtPlusSeqY[NtPos2];

	Codon[0] = g_CharToCompChar[b0];
	Codon[1] = g_CharToCompChar[b1];
	Codon[2] = g_CharToCompChar[b2];

	char aa2 = TranslateCodon3(Codon);
	asserta(aa2 == aa);

	return BaseNr;
	}

void BiPal::Init(const SeqDB &MSA, uint SeqIndexX, uint SeqIndexY)
	{
	m_NtMSA = &MSA;
	m_NtColCount = m_NtMSA->GetColCount();

	m_NtMSA->GetSeqStr(SeqIndexX, m_NtPlusRowX);
	m_NtMSA->GetSeqStr(SeqIndexY, m_NtPlusRowY);

	StripGaps(m_NtPlusRowX, m_NtPlusSeqX);
	StripGaps(m_NtPlusRowY, m_NtPlusSeqY);

	m_LX = SIZE(m_NtPlusSeqX);
	m_LY = SIZE(m_NtPlusSeqY);

	asserta(m_LX%3 == 0);
	asserta(m_LY%3 == 0);

	InitFrames();
	InitOrfs();
	InitNtPlusPosToCol();
	InitAaOrfPosToPosPlus();
	InitAaOrfPosToPosMinus();
	InitAaRowsPlus();
	InitAaRowsMinus();
	ValidateMaps();
	WriteCodons(g_fLog);
	}

void cmd_bipal()
	{
	const string &MSAFN = opt(bipal);

	SeqDB MSA;
	MSA.FromFasta(MSAFN, false);

	BiPal BP;
	BP.Init(MSA, 0, 1);
	}
