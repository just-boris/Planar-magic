object Form1: TForm1
  Left = 127
  Top = 161
  Width = 870
  Height = 640
  Caption = 'POLOZOK'
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object PaintBox1: TPaintBox
    Left = 0
    Top = 0
    Width = 720
    Height = 602
    Align = alClient
  end
  object Panel1: TPanel
    Left = 720
    Top = 0
    Width = 134
    Height = 602
    Align = alRight
    Caption = 'Panel1'
    TabOrder = 0
    object Label1: TLabel
      Left = 25
      Top = 20
      Width = 27
      Height = 13
      Caption = 'Mode'
    end
    object Label2: TLabel
      Left = 66
      Top = 52
      Width = 53
      Height = 13
      Caption = 'N_Substrat'
    end
    object Label3: TLabel
      Left = 72
      Top = 80
      Width = 45
      Height = 13
      Caption = 'Nref_max'
    end
    object Label4: TLabel
      Left = 83
      Top = 112
      Width = 8
      Height = 13
      Caption = 'H'
    end
    object Label5: TLabel
      Left = 80
      Top = 136
      Width = 11
      Height = 13
      Caption = 'W'
    end
    object Label6: TLabel
      Left = 80
      Top = 160
      Width = 20
      Height = 13
      Caption = 'Lam'
    end
    object Axon_3D: TButton
      Left = 32
      Top = 290
      Width = 75
      Height = 25
      Caption = 'Axon_3D'
      TabOrder = 0
      OnClick = Axon_3DClick
    end
    object ScrollBar1: TScrollBar
      Left = 7
      Top = 219
      Width = 120
      Height = 16
      Max = 180
      PageSize = 0
      TabOrder = 1
      OnChange = ScrollBar1Change
    end
    object ScrollBar2: TScrollBar
      Left = 8
      Top = 237
      Width = 121
      Height = 16
      Max = 180
      PageSize = 0
      TabOrder = 2
      OnChange = ScrollBar2Change
    end
    object ScrollBar3: TScrollBar
      Left = 8
      Top = 256
      Width = 121
      Height = 16
      Max = 180
      PageSize = 0
      TabOrder = 3
      OnChange = ScrollBar3Change
    end
    object Strip_waveguide: TButton
      Left = 19
      Top = 183
      Width = 94
      Height = 25
      Caption = 'Strip_waveguide'
      TabOrder = 4
      OnClick = Strip_waveguideClick
    end
    object Edit1: TEdit
      Left = 65
      Top = 16
      Width = 48
      Height = 21
      TabOrder = 5
      Text = 'Edit1'
    end
    object Substrat: TEdit
      Left = 7
      Top = 50
      Width = 54
      Height = 21
      TabOrder = 6
      Text = 'Substrat'
    end
    object Nref_max: TEdit
      Left = 8
      Top = 77
      Width = 54
      Height = 21
      TabOrder = 7
      Text = 'Nref_max'
    end
    object H_waveguide: TEdit
      Left = 9
      Top = 105
      Width = 56
      Height = 21
      TabOrder = 8
      Text = 'H_waveguide'
    end
    object W_waveguide: TEdit
      Left = 8
      Top = 130
      Width = 59
      Height = 21
      TabOrder = 9
      Text = 'W_waveguide'
    end
    object Lam: TEdit
      Left = 10
      Top = 156
      Width = 60
      Height = 21
      TabOrder = 10
      Text = 'Lam'
    end
    object Graph_1d: TButton
      Left = 32
      Top = 320
      Width = 75
      Height = 25
      Caption = 'Graph_1d'
      TabOrder = 11
      OnClick = Graph_1dClick
    end
    object Cylindr: TButton
      Left = 32
      Top = 360
      Width = 75
      Height = 25
      Caption = 'Cylindr'
      TabOrder = 12
      OnClick = CylindrClick
    end
    object Test: TButton
      Left = 24
      Top = 576
      Width = 75
      Height = 25
      Caption = 'Test'
      TabOrder = 13
      OnClick = TestClick
    end
  end
end
