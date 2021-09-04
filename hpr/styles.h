#pragma once
#include "hpr.h"
#include "IPlug_include_in_plug_src.h"

enum Colors
{
	FG1 = 0x005F6B,
	FRAME = 0x68B3AF,
	ON_PRESSED = 0x413E4A,
	BG = 0x3B8183,
	HC	= 0x1B676B,
	HP	= 0x519548,
	HR	= 0x88C425,
	TH	= 0xAACCB1,
	TP	= 0xAACCB1,
	TR	= 0xAACCB1,
	TLC	= 0xFAD089,
	TLP	= 0xFF9C5B,
	TLR	= 0xF5634A,
	INDC = 0xF7E4BE,
};

IColor getColor(Colors c) { return IColor::FromColorCode(c); }

IVStyle baseStyle()
{
	return DEFAULT_STYLE
		.WithColor(kBG, getColor(Colors::BG))
		.WithColor(kHL, COLOR_RED)
		.WithColor(kSH, getColor(Colors::FG1))
		.WithColor(kFR, getColor(Colors::FRAME))
		.WithColor(kX1, getColor(Colors::INDC))
		.WithColor(kX2, COLOR_GREEN)
		.WithColor(kX3, COLOR_BLUE)
		;
}
IVStyle knob1Style()
{
	// Background
	// OFF/Foreground
	// ON/Pressed
	// Frame
	// Highlight
	// Shadow
	// Extra 1
	// Extra 2
	// Extra 3
	IVStyle res = baseStyle()
		.WithColor(kFG, getColor(Colors::HC))
		.WithColor(kPR, getColor(Colors::ON_PRESSED))
		;
	res.labelText = IText(20, getColor(Colors::TLC));
	res.valueText = IText(16, getColor(Colors::TH));

	return res;
}
IVStyle knob2Style()
{
	// Background
	// OFF/Foreground
	// ON/Pressed
	// Frame
	// Highlight
	// Shadow
	// Extra 1
	// Extra 2
	// Extra 3
	IVStyle res = baseStyle()
		.WithColor(kFG, getColor(Colors::HC))
		.WithColor(kPR, getColor(Colors::ON_PRESSED))
		;
	res.labelText = IText(20, getColor(Colors::TLP));
	res.valueText = IText(16, getColor(Colors::TP));

	return res;
}
IVStyle knob3Style()
{
	// Background
	// OFF/Foreground
	// ON/Pressed
	// Frame
	// Highlight
	// Shadow
	// Extra 1
	// Extra 2
	// Extra 3
	IVStyle res = baseStyle()
		.WithColor(kFG, getColor(Colors::HC))
		.WithColor(kPR, getColor(Colors::ON_PRESSED))
		;
	res.labelText = IText(20, getColor(Colors::TLR));
	res.valueText = IText(16, getColor(Colors::TR));

	return res;
}


