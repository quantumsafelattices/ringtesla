/* This is free and unencumbered software released into the public domain.
 *
 * Anyone is free to copy, modify, publish, use, compile, sell, or
 * distribute this software, either in source code form or as a compiled
 * binary, for any purpose, commercial or non-commercial, and by any
 * means.
 *
 * See LICENSE for complete information.
 */

static uint64_t tesla_table[372][3] ={
{0x23D54D594395957,0x6399B2F4FF150616,0x470BE1F464C3A87D},
{0xC891982DA5EBFBE,0x4AF529574B116BC6,0xD501485326F116ED},
{0x7D7C5164E244CB8E,0x52277787D4F693E1,0x1628FC7E9FB8C050},
{0xDA389CE6ACA8222B,0xCDCAC083E38FB7CF,0x1EF736C8514D2CCC},
{0xE101467FFD72E81C,0xEE4059823C2E9CAD,0x27B6938C7215A1A2},
{0x3FD94019551EAC9B,0x01626AC339CF570E,0x3062F85A3B6404AA},
{0x53CE8A56237916B9,0x4C614D535D5D87D0,0x38F865EA521E80BD},
{0xF7B5B7F589C0EEBC,0xEB897EE20BA37C29,0x4172FDA25C995A75},
{0x6BF82DF42C630CC0,0x661A23427A9BF423,0x49CF06CB2B684277},
{0xF27B4E68DBC2B8D6,0x5189E96031A0A739,0x5208F36DD01699A3},
{0x31C0E5185F608DB7,0x722CB0FC95CEE56C,0x5A1D64CD2199CF2B},
{0xD0A4A7D537ADEDE2,0x8F86EE81895423C4,0x62092F72729D263D},
{0xC52E66E1CA81FF26,0x5C98F7541651BDC7,0x69C95EC59C2E6944},
{0x6C4E74E9C734EA92,0x426D5AA1091B0C32,0x715B3829F22EC963},
{0x859895471ED76595,0x1C3C768B076CC1F7,0x78BC3D9B38522840},
{0x6C6DD417C5FBB574,0x2ED427408D2A3D86,0x7FEA2FC7385ECDB3},
{0x0944533846929857,0x3DA87194FF254868,0x86E30FA226B8CD0E},
{0xB5A262F5A6DC406C,0xB47A5E5446417771,0x8DA51F758B1819DA},
{0x3574A6AF61CCF702,0x838B06828369E7DB,0x942EE36AE2D39D1D},
{0x62189D22C7CF7404,0xF3EB041FD6391407,0x9A7F2194A2F9DD0D},
{0xDA41B75704A430FB,0x7E53425E08A3BCFB,0xA094E1799D726C97},
{0x1B27EAB299BC7FCE,0xD7BD8743A77E23C0,0xA66F6B2811070117},
{0x7E030DE9A30D053A,0x938E0A93E8B86168,0xAC0E45D6D47182D7},
{0x35C44FDFF624D41F,0xF858A405E6313BD5,0xB171361C14FB1D2D},
{0xDEC378449B526495,0x101D20C9CD194983,0xB6983BC207132BFF},
{0x8CCE786DA30E817E,0xE31BAF4DCD024C6E,0xBB838F42A875689F},
{0xF16E60D4DFBD3966,0x12C9B172DB1536BF,0xC0339EF44E974B2C},
{0x82D535D36818464B,0xF574B25C97AE8CC5,0xC4A90BF130700B6F},
{0x5FD970951F1A0219,0x79C17233A9136448,0xC8E4A6C4683C4D2D},
{0x53AC4C28E27868F3,0x6C8672CD20E79CFC,0xCCE76BE7113AB6FB},
{0xC0B7FA94AB3C9FB7,0xB778208384D6FB0A,0xD0B2801827C0AAFB},
{0xE06216253C4B93C7,0xB1BF76B96F1A63B8,0xD4472C99B13310D2},
{0xB227F7435B3F240A,0xF5835748B4642A60,0xD7A6DB5D6FAA205A},
{0xE9F81858BF6E9FE8,0xCAA217D26F0EF4BE,0xDAD3132B05016DB6},
{0x1B3B5E4494157332,0xA328C225CCF8824B,0xDDCD73C8EDDB4486},
{0x4EAAE7FC1AC27E10,0x498FF744B4967DC5,0xE097B2312A03BA3A},
{0x889A603D119C478C,0x7B13385DACD9956E,0xE33394D9C02D1E64},
{0x2CA6D32DB8BB8B31,0x3D941A9DC16CE695,0xE5A2F01892DDA2D9},
{0x39AF5B9EA8D4A9CD,0xFD048772254BF007,0xE7E7A2A9374BF087},
{0x5195DDEEF3BA34AD,0xA2ECDA91D31E8D53,0xEA03925AB0880955},
{0xAD1CACAB4B61DEBE,0x56C8982668D48234,0xEBF8A8EA1D43F54F},
{0xD03A503B7EB8FA23,0x7831F41D10291ECA,0xEDC8D10E9053F303},
{0xC913E8DACE344F27,0x516C9C1DDD962987,0xEF75F3B976D53319},
{0xEFE758EF889E383D,0x0DC56A411A96C217,0xF101F58E1DDB19A1},
{0xA18DE8D63E9B9E25,0x644F8C59D20CBAB3,0xF26EB4921F4BFBFB},
{0x680737895D6942C0,0x827FB06EBD534B15,0xF3BE0617BCA60B35},
{0xAF1D1011E0BBE520,0x5A421E31F3175E0C,0xF4F1B4E278D51F93},
{0x8BA8C16DE8F8DBD4,0xA6CB6F27C7114DC1,0xF60B7F8599D382DB},
{0xE39F35AB5CAE4E31,0x6F7A5525B431C3EF,0xF70D16FBA0E9F4DB},
{0x72BD687D734A9B30,0xDC42ED025432BF92,0xF7F81D753E141630},
{0x77679E27BB30054E,0xCE6E92861345ECF3,0xF8CE255DC90D60FA},
{0x5ECC7C3024700983,0xC6E497182DC5FF29,0xF990B092E13C081A},
{0x6F38BF7165800568,0xE2598C1ADE57E4FF,0xFA412FCC7D28C514},
{0x6F780A194EC6AB82,0x41FC435CFBBE7C7E,0xFAE102326A437827},
{0x5132046C72EEAB66,0xC606171F7BC5F609,0xFB71751C05EABFB1},
{0xC16E42FB137E1B9A,0xE79C1442D6163026,0xFBF3C3F6D05DE4BE},
{0x7B1FFEC3ABFE482E,0x582052129160B12E,0xFC6918506091C9F1},
{0xB7C2020E20494C15,0xBC4D2272370850F9,0xFCD28A0033D008C8},
{0x20627AC8EBBB42DC,0x0BD75709FD943130,0xFD311F6DD470942E},
{0x666BF89A705001BC,0x34BF5AB7B19018FE,0xFD85CDEFDEA82E20},
{0xEAF9DBB163F5B583,0xA9A6C3BAE5BABF4A,0xFDD17A4080076AA8},
{0x48A007506F8338C6,0x2D43324512A44AD1,0xFE14F9042C8D7D1D},
{0xC3EFD7A0B39BA5DD,0xE45FDBF656CD46E1,0xFE510F5F69CD1670},
{0x628A6F18F50F9B8B,0x403265E09DD1D2E2,0xFE867398BC429648},
{0x8B2464BCCCC9387D,0xD3109F809D12894D,0xFEB5CDC3F35E8F30},
{0xB833588A7E74144E,0x30619635A8EF9336,0xFEDFB87444BA02F6},
{0xAA6B2C80BD329318,0x8D869FAD26433DC2,0xFF04C172DD4994C8},
{0x36411FE3A7EFA6BF,0x5BD41DD0675DEA43,0xFF256A77C7270829},
{0x602BE4DC34D8388D,0x6174199456958E0D,0xFF4229E33AB9D518},
{0x5E3748CF5550C78E,0xC5E809AA6C0A1BC3,0xFF5B6B75A9CF589E},
{0x0E1D018A15DEBF9E,0x1C3884D88181D1C4,0xFF71910509D2E943},
{0x64EC1CE62D995DDB,0xC72C9017020743F8,0xFF84F32E172A90C1},
{0x8A5A3A059027BB01,0xB64299307CFC518F,0xFF95E2007E3DCE9C},
{0x0EE58CE5BF1572EB,0x57B0F205972755E3,0xFFA4A5A50669F082},
{0xF7E4CC7F59F672FF,0x06F3A212A8F825CF,0xFFB17EFD07CD5C11},
{0x4F9C05963DB48A6D,0x97AFE1A313F73368,0xFFBCA83A9E21EF80},
{0x3FBC1FA258C595EC,0x7D0DACF66B1DCC48,0xFFC6557130AD298A},
{0xC1344EBF83BFC1AE,0x4948272B412909AE,0xFFCEB51E09854E12},
{0x9E8436DE08916FE4,0x9D9C2FE81F3A0420,0xFFD5F0A8D508CFA8},
{0x4C60635571719913,0x17344F013E4C5564,0xFFDC2CDBFD73A742},
{0xBB10743276A70125,0x6E6F17DC50515D29,0xFFE18A54EE14B498},
{0xAC5BD581E2C8FD89,0x658E60E74698D025,0xFFE625EC5DE65F12},
{0x162CB423DAE8046A,0x37738E8295B7D10C,0xFFEA1916D264DF48},
{0xB22525F3B4778F4A,0xB473FB169D5B751A,0xFFED7A3D999C17FB},
{0x56F7143E3EF82E00,0x6C416D1EE02A9F7D,0xFFF05D1085C5F175},
{0x4ACC15D22D33AEEC,0x186F8D6BDB57CC43,0xFFF2D2D0BDAB5B5E},
{0xF3171A26CEDBB7F8,0x53085970E14DFA39,0xFFF4EA94FB89E0F9},
{0xBFDFCE259CD0D024,0xFFB191F42D07E4E3,0xFFF6B18798C313D6},
{0x7ED3C42B10C72AA3,0x7904C91A09CF8EA5,0xFFF8331EC74E5781},
{0x316CA4CE03FE538D,0xB350E41F99E6835E,0xFFF9794F5B0398B7},
{0x60CDC24604B26217,0x6D48116E7B207FDC,0xFFFA8CBA8496A26C},
{0x38AEDE8270067846,0x859483208C2256F2,0xFFFB74D6DEB481A0},
{0xFACE427554154E54,0x68CDD53E36915D04,0xFFFC38152B58B601},
{0xD27B9B1F0772ABC9,0xF2AC52C2149C54E1,0xFFFCDC011C4CEF18},
{0xED5A0A6BF4D27389,0xB4CC8EFEF4F297C9,0xFFFD655E7E05ED58},
{0xFC41DD49270F2524,0xF30D637C903C2A20,0xFFFDD84317D8CA88},
{0x70C6A6703205AD1F,0x1404A5B687CDF10A,0xFFFE382D8FF7497D},
{0x622646EB3F9F19BE,0x98EA23060554FF7D,0xFFFE88199CDDA309},
{0xAA8607B0DC2B37C1,0x5091FAFE768DA296,0xFFFECA91C8FA710F},
{0x8B82CAF2810E5D8D,0x9D4683176DB67755,0xFFFF01BF086FF4DE},
{0x2903B0059A92B9F1,0x09526E856BCAB9A5,0xFFFF2F765BEF9DC9},
{0x3AB1CC001919EE0A,0xFB52A3E1835DC25C,0xFFFF5544B6ECBFE9},
{0x6894C67272A1D3B1,0xF2637858D8FFF0C3,0xFFFF74795AC63875},
{0x8F4CA37BA2E06CAE,0xD06A32EDBBEF4F50,0xFFFF8E2ED419FA3F},
{0x90C1511E889476FD,0x60FF6849EF8DB1DE,0xFFFFA352C3407263},
{0x1FBD1F4740248857,0x2C6EC071720F8680,0xFFFFB4AC94F58CD1},
{0xE1A0BC00D6EC9096,0x2C166A466DD55190,0xFFFFC2E34C7D2F31},
{0xEA3ADE253B071EA0,0xC6A0683B2282A940,0xFFFFCE827D1B07A9},
{0x394A5D06F65B304F,0xF022A7E9EBD9D540,0xFFFFD7FE8D81FA3E},
{0x1C914720455D2CE0,0xD7F665E04E897E73,0xFFFFDFB85CED9F0A},
{0x7A58D6487A54856F,0x6E54D00F84CD7FDA,0xFFFFE6005EE695D2},
{0x75BCC9459361D4F3,0xC37A97D9E22B39D6,0xFFFFEB19403EACAC},
{0x06B27CCA1584F22C,0xD162CFF671F0684A,0xFFFFEF3A25992A7F},
{0xE8BF624FC7FBCC1F,0x57B62D10F1C3B9D6,0xFFFFF29091D2C9F1},
{0xBFD5A8FC0F3B6F21,0x313CEDD46A12EEE3,0xFFFFF54200D0A193},
{0x9ABE682FDF619EAE,0xC76C9B18E26DC01B,0xFFFFF76D41A1BE2E},
{0x25DA6F7CD2D1858A,0x60D5771536D1903C,0xFFFFF92B9970E72C},
{0xBB9B2990D5C6DF8A,0xEF7CA3A2F3FD47B6,0xFFFFFA91B77F38B8},
{0x5D956D614A3EEF17,0xBFD499E76884FF45,0xFFFFFBB081415B08},
{0x02FBDB38EF5F9DF8,0x3CD41439C28F268D,0xFFFFFC95BCBCA78E},
{0x175E7CA4838216CC,0x152A53907980F4C2,0xFFFFFD4C9E643133},
{0x2CF4BD70D4BD7EDA,0x6513CB93C3DF84C2,0xFFFFFDDE3EF43EA5},
{0xCACDA163B1888E52,0x528790415D509FB5,0xFFFFFE51FD22633B},
{0x55DCFB5FBF725022,0x624FE0AFDDCD3033,0xFFFFFEADCE6646B2},
{0x72DDCC919F975639,0x4618D77AA9EA9EC7,0xFFFFFEF681A0CF76},
{0x859D100DF2C94BED,0xAA5F90B60DBBA805,0xFFFFFF2FF5FA75CE},
{0x14671F37DCD09CE4,0xD5E249D501CB00EC,0xFFFFFF5D47F3ED43},
{0x96F57B47111800F5,0x58E6A47D71E01258,0xFFFFFF80F6542514},
{0xE086313FEE6A391D,0xB4F51E537820DFFE,0xFFFFFF9D005A2A1C},
{0xAC1AA427434CE731,0x75A9EE1A030856E9,0xFFFFFFB2FE5F2BA9},
{0xF027CFE04228158A,0x12BEB275CB052937,0xFFFFFFC435E36C06},
{0x8A394B7498229BE4,0xC51502F3428209F9,0xFFFFFFD1A9D707C9},
{0xBD79047BD8100C58,0x6C6BD612353FF8CB,0xFFFFFFDC27CC2B12},
{0xD0B358FCD0EA7D8B,0x593C90A08BBD8B48,0xFFFFFFE452A2928E},
{0x77000CB867FA4578,0x1E92FFEBBA6039F5,0xFFFFFFEAAB234614},
{0x82DA36F3DFFF294C,0xD139001A8372E025,0xFFFFFFEF96EE9D5E},
{0x113DC6D81D5E9B2C,0x3C61770A80E90553,0xFFFFFFF3660D375F},
{0x272239903D82E13D,0xB3E3DE0D39493222,0xFFFFFFF657661106},
{0x09F9A66D9D21D80A,0xF75DFA6AD18376CB,0xFFFFFFF89C4FE423},
{0x534687E69EE0C028,0x128B44B7EB71AFDE,0xFFFFFFFA5B6A0548},
{0x3B8E3B19BD263138,0xCDB846337B08DD87,0xFFFFFFFBB2E0C2C8},
{0xC31282C4E90861ED,0x8CB8B74669BA8E5E,0xFFFFFFFCBA3A85BC},
{0x931D7896F7BAED55,0xAF6C8F4FE42DB8E5,0xFFFFFFFD83C56A37},
{0x060A13760B89FCCC,0x530646DAB6EEC734,0xFFFFFFFE1DB879BD},
{0x7621819382B4F5E5,0x344DB0434026668A,0xFFFFFFFE9317FD96},
{0x89A2AC66B884B301,0xA0A24026A9430382,0xFFFFFFFEEC695985},
{0xF8C764B6CFE2FD22,0xE05D5D80CC667376,0xFFFFFFFF30406B15},
{0x826501CD5C742689,0xF1FC4441685FF827,0xFFFFFFFF63AE6C8E},
{0x243850122819E8AE,0x493A90A0C16D0848,0xFFFFFFFF8A98BF12},
{0x2D421FDA1462B979,0x6A336A5BCE1C8975,0xFFFFFFFFA7FCB370},
{0x1778934730D37496,0x8E17E2904C9E577B,0xFFFFFFFFBE245E55},
{0x3793F09D55C72F17,0x0E1B7A911AA4BEA4,0xFFFFFFFFCECFAE35},
{0x9B55519C9E14F716,0x205D4A9DDDD0E428,0xFFFFFFFFDB544DE7},
{0xAB065107CCAAEE5B,0xC74D48EF86350B2E,0xFFFFFFFFE4B65606},
{0x97472C88D79C56AC,0xB2D8E3499053EAE4,0xFFFFFFFFEBBB6213},
{0xB0D9E4FF0EA7EF06,0xCE5737655B097765,0xFFFFFFFFF0F947B7},
{0x1908DA40619B0E55,0xD59407530A06A6C9,0xFFFFFFFFF4E169E2},
{0x7F28AABAC1E1D786,0x7F75BAA512E42727,0xFFFFFFFFF7C96B1F},
{0x298DCD7CADF23C39,0x20687B0B5B8602D9,0xFFFFFFFFF9F1D7AE},
{0xCDD7AACE5CCFADEF,0xFCA3890B96F795E0,0xFFFFFFFFFB8B3E38},
{0x49A340C95F4F6FD6,0x5D37C6D1028C7921,0xFFFFFFFFFCBA137D},
{0x7C4028961DF096DE,0xEF345AC3D286D640,0xFFFFFFFFFD99A887},
{0xED8AF462793EF461,0x83203DC5075A4E43,0xFFFFFFFFFE3E6AEE},
{0xB82281600EAFAB9E,0xF86DDE4ACF2EC3F0,0xFFFFFFFFFEB799CE},
{0xFC07D8FFECB575B2,0x68CA249105399584,0xFFFFFFFFFF109087},
{0x781047D91AC0D66B,0x653E4202318FA7DD,0xFFFFFFFFFF51C07B},
{0xFDB906A7E42FFF6B,0x96FC0999FD8039F4,0xFFFFFFFFFF816D52},
{0x52E39E19D6ACEA13,0x866D409749E4B63C,0xFFFFFFFFFFA43A8D},
{0x07F32ECEE18A326C,0x7AB88C045043D742,0xFFFFFFFFFFBD95C9},
{0x2B2CD71BB8430A10,0x81D589BA0D0F843B,0xFFFFFFFFFFD00658},
{0xC5F8A97372DA3ECE,0xE7301FCB612AB4D4,0xFFFFFFFFFFDD68BF},
{0xCEF039F06FF0800C,0x176DDE0D85AC1FCA,0xFFFFFFFFFFE71B19},
{0xF08CE10A3C8ACCC7,0xAA5D0FFBADBC118E,0xFFFFFFFFFFEE1E1A},
{0x5CB6594FCC3F3CBE,0xE1609E1F240A79B3,0xFFFFFFFFFFF32DA2},
{0x9C0FC18FAF846554,0x2C27E2BE6234EE35,0xFFFFFFFFFFF6D2E9},
{0x0B39BFF3D26DE554,0xEEE83B91544EC8FA,0xFFFFFFFFFFF971F3},
{0xBB2DFD915CAEDBC4,0x95589BF05B6BABFA,0xFFFFFFFFFFFB5388},
{0x26C7E425384AA85A,0x1173CAFCA095A177,0xFFFFFFFFFFFCAC7E},
{0xE103B12F610543FD,0x2E7D105EFA647BCD,0xFFFFFFFFFFFDA31F},
{0x1B7F8F9FF38CBBD4,0xFCC045FC427F790A,0xFFFFFFFFFFFE531D},
{0x360AA5F5C99F61AA,0x7DF48C56171E53DC,0xFFFFFFFFFFFED078},
{0xABFF4DF715522CCA,0xF19D3C638C6DB252,0xFFFFFFFFFFFF2995},
{0x6B85CE4CABDDD8D7,0x9E7882EE506D9B2E,0xFFFFFFFFFFFF68D1},
{0x7F24BFA52FD52327,0x2B1FBF64F6C95D13,0xFFFFFFFFFFFF959A},
{0xF7FCCD0727B29026,0x5691E291C54E9729,0xFFFFFFFFFFFFB542},
{0x05A8440EFE26FFEE,0x458370A1989D45C4,0xFFFFFFFFFFFFCB98},
{0xA59472BFDB874433,0xF74E2EC102CD9563,0xFFFFFFFFFFFFDB52},
{0x30F4D598B1E5C22A,0x4C293E22717D46C2,0xFFFFFFFFFFFFE661},
{0x64B9E5DB9B551E86,0xF530826AD4A829B7,0xFFFFFFFFFFFFEE22},
{0x6EAC72B364AF9844,0x654F50E393C69250,0xFFFFFFFFFFFFF391},
{0x4AE596CC4DF75470,0x325A6FA33AA9E6D3,0xFFFFFFFFFFFFF75D},
{0xF769B2CF0CD0738D,0x1F4FF57C13B3FA43,0xFFFFFFFFFFFFFA03},
{0x8BEB16815B83875C,0x25EF863DFF5C4975,0xFFFFFFFFFFFFFBDB},
{0x2F4ECE4B47528438,0x30853DA5DA72EE10,0xFFFFFFFFFFFFFD23},
{0x125570CB90234707,0xBC74CF4585941EB7,0xFFFFFFFFFFFFFE06},
{0x00F2272222C19F0E,0x46BB9CA0332E2998,0xFFFFFFFFFFFFFEA4},
{0x70FBCEC1A2AA20C0,0x245AA2921CCDAAEE,0xFFFFFFFFFFFFFF11},
{0x65AF3CFE5690CA33,0x3ADBA46CEA9DE851,0xFFFFFFFFFFFFFF5C},
{0x0A6F10EE29D05B3F,0xEC26435CEADAE5AC,0xFFFFFFFFFFFFFF8F},
{0xE96A9422E7C5085A,0x711996BCC7393656,0xFFFFFFFFFFFFFFB3},
{0xB5D414503A1B18B1,0xCD4F9B72DCE976AD,0xFFFFFFFFFFFFFFCB},
{0x4523FEC3CC5B7121,0x7A358FAB598D1F37,0xFFFFFFFFFFFFFFDC},
{0x1849644B62B099DF,0xDEE57291616E155B,0xFFFFFFFFFFFFFFE7},
{0x7AAEF8E77B0B5008,0xA3EEF45B28496786,0xFFFFFFFFFFFFFFEF},
{0xF0A36ED6D4F16501,0xEDCC1679F8C16C47,0xFFFFFFFFFFFFFFF4},
{0xF916E58C90F004C9,0x8594EBCE15BFD8AD,0xFFFFFFFFFFFFFFF8},
{0xED35DF335CC6D5A4,0xF5484CA3CA2C2FA5,0xFFFFFFFFFFFFFFFA},
{0xD012960DF21E4015,0x9B695FFC5D26056A,0xFFFFFFFFFFFFFFFC},
{0x8417958E19BFD6CB,0xB8932C4093C7785C,0xFFFFFFFFFFFFFFFD},
{0x4E8684CA58A50B46,0x78D9586BCC11CD62,0xFFFFFFFFFFFFFFFE},
{0xE2B0188EEC4F2A3B,0xFA3F31F6869230E3,0xFFFFFFFFFFFFFFFE},
{0x825440FC484E96D9,0x512A564FD67BE67A,0xFFFFFFFFFFFFFFFF},
{0x7E346D0E09CE1BD2,0x8B708B22011523DB,0xFFFFFFFFFFFFFFFF},
{0x101A76E693CFF0CF,0xB26F847C3F6CD018,0xFFFFFFFFFFFFFFFF},
{0x9E9AA26DCE6C9AD4,0xCC7B503F8982722A,0xFFFFFFFFFFFFFFFF},
{0x59B3EB8DA1749C38,0xDDD86F7CB3F159A3,0xFFFFFFFFFFFFFFFF},
{0x16741CA9D4080934,0xE966341689E5D9FC,0xFFFFFFFFFFFFFFFF},
{0x8228CA6F2BA83D7F,0xF11293192EFD5166,0xFFFFFFFFFFFFFFFF},
{0xA2EE123887417A96,0xF628BA4399A8D9EE,0xFFFFFFFFFFFFFFFF},
{0x2130B24DA844624B,0xF9864691566151BF,0xFFFFFFFFFFFFFFFF},
{0xE4CA88DDC32EB7AA,0xFBBF3A6A633FB8AC,0xFFFFFFFFFFFFFFFF},
{0x8754CF86468E8C7F,0xFD363EE82246CB24,0xFFFFFFFFFFFFFFFF},
{0x8E56FEF8707416C6,0xFE2CF77D599C5C7F,0xFFFFFFFFFFFFFFFF},
{0x44041B97A1E0BDCF,0xFECEF9DA8CEB417E,0xFFFFFFFFFFFFFFFF},
{0xD998055BF98E113C,0xFF39288A2ED2012F,0xFFFFFFFFFFFFFFFF},
{0x5F6999A3EB09443E,0xFF7E9EAE14C5D553,0xFFFFFFFFFFFFFFFF},
{0x7AD47611C5978CFC,0xFFABF94810ADC0C7,0xFFFFFFFFFFFFFFFF},
{0xE175A3117AE1C970,0xFFC987F1FCCA3F89,0xFFFFFFFFFFFFFFFF},
{0xD7563BE36C09DF58,0xFFDCC1D9876AB158,0xFFFFFFFFFFFFFFFF},
{0x00C898360113902E,0xFFE93D5D89BAE166,0xFFFFFFFFFFFFFFFF},
{0xC8DDA434D6A3A289,0xFFF1541DA252D2BC,0xFFFFFFFFFFFFFFFF},
{0x55B97648FCFF1E12,0xFFF68F76D5621AEE,0xFFFFFFFFFFFFFFFF},
{0x49F2ED5AF04AB4CE,0xFFF9F01C418F3899,0xFFFFFFFFFFFFFFFF},
{0xE31F5C1482A0C773,0xFFFC1D3F6DD93BB9,0xFFFFFFFFFFFFFFFF},
{0x1327F5666CA55398,0xFFFD83901E02FEF2,0xFFFFFFFFFFFFFFFF},
{0x80BE62C979788293,0xFFFE6992BC776A35,0xFFFFFFFFFFFFFFFF},
{0x8E41EBF252585B07,0xFFFEFCF15F456189,0xFFFFFFFFFFFFFFFF},
{0x8D50149381064402,0xFFFF5B2F865BA43E,0xFFFFFFFFFFFFFFFF},
{0x2CF41FE4041144DD,0xFFFF97570A9F823A,0xFFFFFFFFFFFFFFFF},
{0x332D4D06CA213969,0xFFFFBDA9DB8C01E5,0xFFFFFFFFFFFFFFFF},
{0x6B276480AEB6A0BB,0xFFFFD60863022046,0xFFFFFFFFFFFFFFFF},
{0x683919B83D924679,0xFFFFE57FDE4022B9,0xFFFFFFFFFFFFFFFF},
{0x3378A075A6FCBB4A,0xFFFFEF4C221E3073,0xFFFFFFFFFFFFFFFF},
{0x45097500F2407905,0xFFFFF57E138AD48E,0xFFFFFFFFFFFFFFFF},
{0xEA44300AB1DBBF33,0xFFFFF966F417C302,0xFFFFFFFFFFFFFFFF},
{0x94A42A620BF2E767,0xFFFFFBDD68292CA2,0xFFFFFFFFFFFFFFFF},
{0x9008C635ED05B33F,0xFFFFFD69C7AF981A,0xFFFFFFFFFFFFFFFF},
{0x24F1650A790868F8,0xFFFFFE62837DFFC2,0xFFFFFFFFFFFFFFFF},
{0xD0F01E4D3DF0E3CF,0xFFFFFEFE4E1E3D32,0xFFFFFFFFFFFFFFFF},
{0x28AC57D6D9E3F487,0xFFFFFF5FB2F9E906,0xFFFFFFFFFFFFFFFF},
{0xD62AE6BDFC754A20,0xFFFFFF9C787644C7,0xFFFFFFFFFFFFFFFF},
{0xA89E1EFBAB379C61,0xFFFFFFC2519D491C,0xFFFFFFFFFFFFFFFF},
{0x4EF121CC3355A344,0xFFFFFFD9D8975816,0xFFFFFFFFFFFFFFFF},
{0x98D21DDA717F5655,0xFFFFFFE87186DED9,0xFFFFFFFFFFFFFFFF},
{0xC7D7ACC6FF9CB298,0xFFFFFFF17BB79904,0xFFFFFFFFFFFFFFFF},
{0x5439A93F576D704A,0xFFFFFFF712216BB5,0xFFFFFFFFFFFFFFFF},
{0x02F0AAD8F43A851A,0xFFFFFFFA849B5D80,0xFFFFFFFFFFFFFFFF},
{0xCD1D7F00015CA311,0xFFFFFFFCA4030A25,0xFFFFFFFFFFFFFFFF},
{0x2D2C452B1090A96C,0xFFFFFFFDF1FE7479,0xFFFFFFFFFFFFFFFF},
{0x4112EF74F8770152,0xFFFFFFFEBEE0006F,0xFFFFFFFFFFFFFFFF},
{0xB8462ADF1139FD7E,0xFFFFFFFF3C525E4A,0xFFFFFFFFFFFFFFFF},
{0x9E51B6165D25D51C,0xFFFFFFFF88FC8F10,0xFFFFFFFFFFFFFFFF},
{0x0C692D8B9556B673,0xFFFFFFFFB7C02510,0xFFFFFFFFFFFFFFFF},
{0x6307A4CFE196A7E5,0xFFFFFFFFD438C86F,0xFFFFFFFFFFFFFFFF},
{0x339DD46424A3632A,0xFFFFFFFFE585EE24,0xFFFFFFFFFFFFFFFF},
{0xB8B85AC01222FD1D,0xFFFFFFFFF0045807,0xFFFFFFFFFFFFFFFF},
{0xC0B1051D7C95771D,0xFFFFFFFFF65EA8AB,0xFFFFFFFFFFFFFFFF},
{0x418B1BC9FA3CE686,0xFFFFFFFFFA3558B3,0xFFFFFFFFFFFFFFFF},
{0xD149B722B60EAADB,0xFFFFFFFFFC86018E,0xFFFFFFFFFFFFFFFF},
{0x83E46DAE7989EA13,0xFFFFFFFFFDEAC3B4,0xFFFFFFFFFFFFFFFF},
{0xDCFFBA2C4BB7E48F,0xFFFFFFFFFEC11D20,0xFFFFFFFFFFFFFFFF},
{0x2ED0966DAD63CD5C,0xFFFFFFFFFF41A82E,0xFFFFFFFFFFFFFFFF},
{0x51F9935C124819DF,0xFFFFFFFFFF8E98F2,0xFFFFFFFFFFFFFFFF},
{0x7F7C3E088734C277,0xFFFFFFFFFFBC9054,0xFFFFFFFFFFFFFFFF},
{0x376697885DBA068F,0xFFFFFFFFFFD7F935,0xFFFFFFFFFFFFFFFF},
{0x585E19F91C0F3431,0xFFFFFFFFFFE84976,0xFFFFFFFFFFFFFFFF},
{0x60CD7C69FE529ABD,0xFFFFFFFFFFF1FA5C,0xFFFFFFFFFFFFFFFF},
{0x24F06245AC0EE7B5,0xFFFFFFFFFFF7B95A,0xFFFFFFFFFFFFFFFF},
{0x24F693CBD275F878,0xFFFFFFFFFFFB1FE5,0xFFFFFFFFFFFFFFFF},
{0x954476B69EA576CD,0xFFFFFFFFFFFD221D,0xFFFFFFFFFFFFFFFF},
{0x55FA7375F4005498,0xFFFFFFFFFFFE5149,0xFFFFFFFFFFFFFFFF},
{0xC635D0E8A897B657,0xFFFFFFFFFFFF03B0,0xFFFFFFFFFFFFFFFF},
{0xD8F19730195273EE,0xFFFFFFFFFFFF6C79,0xFFFFFFFFFFFFFFFF},
{0xC573393430A0A4A8,0xFFFFFFFFFFFFA9E7,0xFFFFFFFFFFFFFFFF},
{0x858DB21F585DFE7A,0xFFFFFFFFFFFFCDD9,0xFFFFFFFFFFFFFFFF},
{0x98DC43C1B852D461,0xFFFFFFFFFFFFE2D7,0xFFFFFFFFFFFFFFFF},
{0x470C1393F2184929,0xFFFFFFFFFFFFEF14,0xFFFFFFFFFFFFFFFF},
{0xF1EED20F08A29055,0xFFFFFFFFFFFFF632,0xFFFFFFFFFFFFFFFF},
{0x69CFA475D53882E4,0xFFFFFFFFFFFFFA55,0xFFFFFFFFFFFFFFFF},
{0xEE6708E3B78D07DB,0xFFFFFFFFFFFFFCBA,0xFFFFFFFFFFFFFFFF},
{0xDF33C1323CB1E429,0xFFFFFFFFFFFFFE1D,0xFFFFFFFFFFFFFFFF},
{0xD43EDBB079FD4F0A,0xFFFFFFFFFFFFFEEA,0xFFFFFFFFFFFFFFFF},
{0xF4CAF4656275C2DF,0xFFFFFFFFFFFFFF60,0xFFFFFFFFFFFFFFFF},
{0xE9062B8AD6AE6B21,0xFFFFFFFFFFFFFFA4,0xFFFFFFFFFFFFFFFF},
{0xED82F58A28572373,0xFFFFFFFFFFFFFFCB,0xFFFFFFFFFFFFFFFF},
{0x49D2E11F49E6DADA,0xFFFFFFFFFFFFFFE2,0xFFFFFFFFFFFFFFFF},
{0x142E10532D8906EE,0xFFFFFFFFFFFFFFEF,0xFFFFFFFFFFFFFFFF},
{0x619A53688ADE3BF7,0xFFFFFFFFFFFFFFF6,0xFFFFFFFFFFFFFFFF},
{0x8AE4CCBB76F7AE35,0xFFFFFFFFFFFFFFFA,0xFFFFFFFFFFFFFFFF},
{0xE8CDAFB9B843D803,0xFFFFFFFFFFFFFFFC,0xFFFFFFFFFFFFFFFF},
{0x40C73C9EAE449D16,0xFFFFFFFFFFFFFFFE,0xFFFFFFFFFFFFFFFF},
{0x03AF13E8F42A12E7,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0x71EA0B9FAF02B6BA,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xB0235A99BE80DB67,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xD3324C4AE001CE02,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xE6E970E66423EC72,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xF1FA6367DEFD0E08,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xF82D82B9BC8677FE,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFBA4FF0DFB2E306D,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFD945290A43C78DF,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFEA840DEC3F1663F,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFF41ACF6DC7A75E0,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFF96D1E517169B58,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFC5FB9FAB71B846,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFE00EDBADDD75B2,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFEE727158615EA6,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFF65F4F936A4712,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFABAAC593EA300,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFD1EAE6FDCFD9D,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFE6DDC63BDE716,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFF25151FE1560D,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFF890C946980A9,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFBF7D08CBE8B7,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFDD1423ABFE53,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFED21C820D6C6,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFF5D321E0BE06,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFFA85DF66E0CA,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFFD0EB15CB106,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFFE6C09F29DDD,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFFF27C8711247,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFFF8C7D3599E1,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFFFC267EBC9D7,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFFFDF37C51E18,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFFFEE95CF0CF6,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFFFF6C422840B,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFFFFB1CF71797,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFFFFD6B28059A,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFFFFEA392E9CB,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFFFFF48A3C3F4,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFFFFF9FB002C6,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFFFFFCD8176AD,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFFFFFE593ACCC,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFFFFFF232F43D,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFFFFFF8CE23D8,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFFFFFFC4199D1,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFFFFFFE0E3DD2,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFFFFFFEFDF7D4,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFFFFFFF7A7DB1,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFFFFFFFBB0CF7,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFFFFFFFDC7574,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFFFFFFFEDB6CF,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFFFFFFFF69C10,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFFFFFFFFB2FD3,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFFFFFFFFD899E,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFFFFFFFFEBE1A,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFFFFFFFFF5BEF,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFFFFFFFFFAC89,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFFFFFFFFFD59E,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFFFFFFFFFEA85,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFFFFFFFFFF522,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFFFFFFFFFFA83,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFFFFFFFFFFD3C,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFFFFFFFFFFE9C,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFFFFFFFFFFF4D,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFFFFFFFFFFFA6,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFFFFFFFFFFFD3,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFFFFFFFFFFFE9,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFFFFFFFFFFFF4,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFFFFFFFFFFFFA,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFFFFFFFFFFFFD,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFFFFFFFFFFFFE,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF},
{0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF}
} ;
