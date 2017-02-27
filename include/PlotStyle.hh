class TStyle;

class PlotStyle {
public:
  static void setTklayoutStyle();
  static void setTDRStyle();
private:
  static TStyle *tdrStyle;
  static TStyle *tkLayoutStyle;
  void tdrGrid(bool gridOn);
  void fixOverlay();
};
