//---------------------------------------------------------------------------

#ifndef Poloz1H
#define Poloz1H
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <ExtCtrls.hpp>
//---------------------------------------------------------------------------
class TForm1 : public TForm
{
__published:	// IDE-managed Components
        TPanel *Panel1;
        TPaintBox *PaintBox1;
        TButton *Axon_3D;
        TScrollBar *ScrollBar1;
        TScrollBar *ScrollBar2;
        TScrollBar *ScrollBar3;
        TButton *Strip_waveguide;
        TEdit *Edit1;
        TLabel *Label1;
        TEdit *Substrat;
        TEdit *Nref_max;
        TEdit *H_waveguide;
        TEdit *W_waveguide;
        TEdit *Lam;
        TLabel *Label2;
        TLabel *Label3;
        TLabel *Label4;
        TLabel *Label5;
        TLabel *Label6;
        TButton *Graph_1d;
        TButton *Cylindr;
        TButton *Test;
        void __fastcall Axon_3DClick(TObject *Sender);
        void __fastcall ScrollBar1Change(TObject *Sender);
        void __fastcall ScrollBar2Change(TObject *Sender);
        void __fastcall ScrollBar3Change(TObject *Sender);
        void __fastcall Strip_waveguideClick(TObject *Sender);
        void __fastcall Graph_1dClick(TObject *Sender);
        void __fastcall CylindrClick(TObject *Sender);
        void __fastcall TestClick(TObject *Sender);
private:	// User declarations
public:		// User declarations
        __fastcall TForm1(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TForm1 *Form1;
//---------------------------------------------------------------------------
#endif
