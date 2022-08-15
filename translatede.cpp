#include "myutils.h"
#include "fastaseqsource.h"
#include "seqinfo.h"
#include "objmgr.h"
#include "translator.h"
#include "alpha.h"

static void LogRows(const string &NtRow, const string &AaRow, uint W)
	{
	uint ColCount = SIZE(NtRow);
	asserta(SIZE(AaRow) == ColCount);

//	const uint W = 4*20;
	uint StartCol = 0;
	for (;;)
		{
		uint EndCol = StartCol + W - 1;
		if (EndCol >= ColCount)
			EndCol = ColCount - 1;

		for (uint Col = StartCol; Col <= EndCol; ++Col)
			Log("%c", NtRow[Col]);
		Log("\n");

		for (uint Col = StartCol; Col <= EndCol; ++Col)
			Log("%c", AaRow[Col]);
		Log("\n");

		StartCol = EndCol + 1;
		if (EndCol >= ColCount - 1)
			break;
		Log("\n");
		}
	}


void TranslateDE(const byte *Seq, uint L)
	{
	asserta(L%3 == 0);
	string NtRow;
	string AaRow;
	for (uint Start = 0; Start + 5 < L; Start += 3)
		{
		string Wmer;
		for (uint i = 0; i < 5; ++i)
			Wmer.push_back(Seq[Start+i]);

		NtRow += ' ';
		NtRow += Wmer;

		char PlusCodon[3];
		char MinusCodon[3];

		PlusCodon[0] = Wmer[0];
		PlusCodon[1] = Wmer[1];
		PlusCodon[2] = Wmer[2];

		MinusCodon[0] = g_CharToCompChar[Wmer[4]];
		MinusCodon[1] = g_CharToCompChar[Wmer[3]];
		MinusCodon[2] = g_CharToCompChar[Wmer[2]];

		char AaPlus = TranslateCodon3(PlusCodon);
		char AaMinus = TranslateCodon3(MinusCodon);

		AaRow += "  ";
		AaRow += AaPlus;
		AaRow += " ";
		AaRow += AaMinus;
		AaRow += " ";
		}

	Log("\n");
	//Log("%s\n", NtRow.c_str());
	//Log("%s\n", AaRow.c_str());
	LogRows(NtRow, AaRow, (80/6)*6);
	}

void TranslateDE(const string &Seq)
	{
	TranslateDE((const byte *) Seq.c_str(), SIZE(Seq));
	}

void LogDE(const char *Label, const byte *Seq, uint L)
	{
	Translator T;

	T.SetSeq(Seq, L, true);
	int PlusFrame = 0;
	int MinusFrame = 0;
	T.GetEndlessFrames2(PlusFrame, MinusFrame);
	uint Phase = GetPhase(PlusFrame, MinusFrame);
	if (Phase != 3)
		{
		Log(">%s\n", Label);
		Log("(DE phase %u)\n", Phase);
		return;
		}
	if (!(PlusFrame > 0 && MinusFrame < 0))
		{
		Log(">%s\n", Label);
		Log("(Not DE)\n");
		return;
		}

	string AaPlus;
	string NtPlus;
	string AaMinus;
	string NtMinus;

	T.SetFrame(PlusFrame);
	T.GetAaStr(AaPlus);
	T.GetNtStr(NtPlus);
		
	T.SetFrame(MinusFrame);
	T.GetAaStr(AaMinus);
	T.GetNtStr(NtMinus);

	string AaPlusRow;
	string NtPlusRow;
	string AaMinusRow;
	string NtMinusRow;

	for (uint i = 0; i < SIZE(AaPlus); ++i)
		{
		NtPlusRow += ' ';
		NtPlusRow += NtPlus[3*i];
		NtPlusRow += NtPlus[3*i+1];
		NtPlusRow += NtPlus[3*i+2];

		AaPlusRow += "   ";
		AaPlusRow += AaPlus[i];
		}

	for (uint i = 0; i < SIZE(AaMinus); ++i)
		{
		NtMinusRow += ' ';
		NtMinusRow += NtMinus[3*i];
		NtMinusRow += NtMinus[3*i+1];
		NtMinusRow += NtMinus[3*i+2];

		AaMinusRow += "   ";
		AaMinusRow += AaMinus[i];
		}

	Log("\n_______________________________________________________\n");
	Log(">%s\n", Label);
	Log("\n");
	Log("Plus:\n");
	LogRows(NtPlusRow, AaPlusRow, (80/4)*4);
	Log("\n");
	Log("Minus:\n");
	LogRows(NtMinusRow, AaMinusRow, (80/4)*4);

	TranslateDE(NtPlus);
	}

void cmd_translatede()
	{
	FASTASeqSource SS;
	SS.Open(opt(translatede));

	const uint DictSize = 1024; // myipow(4, 5);
	vector<uint> Counts(DictSize, 0);

	SeqInfo *SI = ObjMgr::GetSeqInfo();
	for (;;)
		{
		bool Ok = SS.GetNext(SI);
		if (!Ok)
			break;

		LogDE(SI->m_Label, SI->m_Seq, SI->m_L);
		}
	}
