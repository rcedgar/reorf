#pragma once

#include "seqdb.h"
class Translator;

class BiPal
	{
public:
	Translator m_Translator;

	const SeqDB *m_NtMSA = 0;

	uint m_SeqIndexX = UINT_MAX;
	uint m_SeqIndexY = UINT_MAX;

	uint m_NtColCount = UINT_MAX;

	uint m_LX = UINT_MAX;
	uint m_LY = UINT_MAX;

	int m_PlusFrame = 0;
	int m_MinusFrame = 0;

	uint m_PlusFrameOffset = UINT_MAX;
	uint m_MinusFrameOffset = UINT_MAX;

// Aligned nt sequences from m_NtMSA
	string m_NtPlusRowX;
	string m_NtPlusRowY;

// Unaligned nt sequences from m_NtMSA
	string m_NtPlusSeqX;
	string m_NtPlusSeqY;

// Plus-strand unaligned nt sequences
// n=(m_PlusFrame-1) letters deleted from start, added to end
// compared to nt sequence from m_NtMSA
	string m_NtPlusOrfX;
	string m_NtPlusOrfY;

	string m_NtMinusOrfX;
	string m_NtMinusOrfY;

// Translated m_NtPlusOrfX/Y
	string m_AaPlusOrfX;
	string m_AaPlusOrfY;

	string m_AaMinusOrfX;
	string m_AaMinusOrfY;

	vector<uint> m_AaPlusOrfPosXToPosY;
	vector<uint> m_AaMinusOrfPosXToPosY;

// Aligned aa sequences
	string m_AaPlusRowX;
	string m_AaPlusRowY;

	string m_AaMinusRowX;
	string m_AaMinusRowY;

// PosToCol mappings
	vector<uint> m_NtPlusPosToColX;
	vector<uint> m_NtPlusPosToColY;

	vector<uint> m_AaPlusPosToColX;
	vector<uint> m_AaPlusPosToColY;

	vector<uint> m_AaMinusPosToColX;
	vector<uint> m_AaMinusPosToColY;

// ColToPos mappings
	vector<uint> m_NtPlusColToPosX;
	vector<uint> m_NtPlusColToPosY;

	vector<uint> m_AaPlusColToPosX;
	vector<uint> m_AaPlusColToPosY;

	vector<uint> m_AaMinusColToPosX;
	vector<uint> m_AaMinusColToPosY;

public:
	void Init(const SeqDB &MSA, uint SeqIndexX, uint SeqIndexY);

	uint CalcAaOrfLengthX() const;
	uint CalcAaOrfLengthY() const;

	uint PlusOrfAaPosToPlusNtPosX(uint PlusOrfAaPos, uint BaseNr) const;
	uint PlusOrfAaPosToPlusNtPosY(uint PlusOrfAaPos, uint BaseNr) const;

	uint MinusOrfAaPosToPlusNtPosX(uint MinusOrfAaPos, uint BaseNr) const;
	uint MinusOrfAaPosToPlusNtPosY(uint MinusOrfAaPos, uint BaseNr) const;

	uint PlusNtPosToPlusOrfAaPosX(uint PlusNtPosX, uint &BaseNr) const;
	uint PlusNtPosToPlusOrfAaPosY(uint PlusNtPosY, uint &BaseNr) const;

	uint PlusNtPosToMinusOrfAaPosX(uint PlusNtPosX, uint &BaseNr) const;
	uint PlusNtPosToMinusOrfAaPosY(uint PlusNtPosY, uint &BaseNr) const;

	void LogMaps() const;
	void WriteCodons(FILE *f) const;

	void ValidateMaps() const;

	uint GetCodonPlusX(uint PlusNtPosX, char Codon[3]) const;
	uint GetCodonPlusY(uint PlusNtPosY, char Codon[3]) const;

	uint GetCodonMinusX(uint PlusNtPosX, char Codon[3]) const;
	uint GetCodonMinusY(uint PlusNtPosY, char Codon[3]) const;

	void GetCodonInfoPlusMinusX(
	  uint PlusNtPosX, char *PlusSeq,
	  char &AaPlus, char &AaMinus,
	  uint &Offset, uint &StartPlus, uint &StartMinus) const;

	void GetCodonInfoPlusMinusY(
	  uint PlusNtPosY, char *PlusSeq,
	  char &AaPlus, char &AaMinus,
	  uint &Offset, uint &StartPlus, uint &StartMinus) const;

	void WriteCodonInfo(FILE *f, const char *PlusSeq, uint Offset,
	  uint StartPlus, uint StartMinus) const;

	void GetBarStr(const char *PlusSeq, uint Offset, string &BarStr) const;
	void GetCodonPlusStr(const char *PlusSeq, uint Start, string &Str) const;
	void GetCodonMinusStr(const char *PlusSeq, uint Start, string &Str) const;

private:
	void InitFrames();
	void InitOrfs();
	void InitNtPlusPosToCol();
	void InitAaOrfPosToPosPlus();
	void InitAaOrfPosToPosMinus();
	void InitAaRowsPlus();
	void InitAaRowsMinus();
	};
