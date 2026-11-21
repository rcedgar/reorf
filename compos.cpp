#include "myutils.h"
#include "fastaseqsource.h"
#include "objmgr.h"
#include "seqinfo.h"
#include <ctype.h>

static void TruncLabel(const char *Defline, string &Label)
	{
	Label.clear();
	if (Defline == 0)
		return;
	for (const char *p = Defline; *p; ++p)
		{
		unsigned char c = (unsigned char) *p;
		if (isspace(c))
			break;
		Label.push_back((char) c);
		}
	}

void cmd_compos()
	{
	if (g_Arg1.empty())
		Die("Missing FASTA file");
	if (!optset_output || opt(output).empty())
		Die("Missing -output");

	FASTASeqSource SS;
	SS.Open(g_Arg1);

	FILE *f = CreateStdioFile(opt(output));

	SeqInfo *SI = ObjMgr::GetSeqInfo();
	for (;;)
		{
		bool Ok = SS.GetNext(SI);
		if (!Ok)
			break;

		uint NA = 0;
		uint NC = 0;
		uint NG = 0;
		uint NT = 0;
		const byte *Seq = SI->m_Seq;
		const uint L = SI->m_L;
		for (uint i = 0; i < L; ++i)
			{
			byte c = (byte) mytoupper(Seq[i]);
			if (c == 'A')
				++NA;
			else if (c == 'C')
				++NC;
			else if (c == 'G')
				++NG;
			else if (c == 'T')
				++NT;
			}

		string Label;
		TruncLabel(SI->m_Label, Label);
		Pf(f, "%s\t%u\t%u\t%u\t%u\t%u\n",
		  Label.c_str(), L, NA, NC, NG, NT);
		}

	SS.Close();
	CloseStdioFile(f);
	}
