import sys
from PyQt5.QtWidgets import (
    QApplication, QWidget, QLabel, QLineEdit, QPushButton, QTableView,
            QVBoxLayout, QHBoxLayout, QFrame, QStackedWidget, QTextEdit, QDialog, QScrollArea, QCheckBox, QMessageBox, QSizePolicy
)
from PyQt5.QtGui import QStandardItemModel, QStandardItem
from PyQt5.QtCore import Qt, QSortFilterProxyModel
from PyQt5.QtGui import QPixmap, QFont
from io import BytesIO
from PyQt5.QtCore import Qt
from PyQt5.QtWebEngineWidgets import QWebEngineView

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas, NavigationToolbar2QT
import matplotlib.pyplot as plt
import inspect
import py3Dmol

### modules ###
#import protein_analysis   # practicemodule
from backend import g1_protein
from backend import g3_variant



# -----------------------------------------------------------
# PAGE 1: WELCOME PAGE
# -----------------------------------------------------------
class WelcomePage(QWidget):
    def __init__(self, parent):
        super().__init__()
        self.parent = parent  #cce7ff
        self.setStyleSheet("""
            QWidget {
                background: qlineargradient(
                    x1:0, y1:0, x2:1, y2:1,
                    stop:0 #d9ecff,
                    stop:1 #f4f9ff
                );
            }
        """)
        ### Title ###
        title = QLabel("ProSeqNet")
        title.setFont(QFont("Segoe UI", 24, QFont.Bold))
        title.setAlignment(Qt.AlignCenter)

        
        subtitle = QLabel("Protein Sequence, Variant & Network Analysis")
        subtitle.setFont(QFont("Segoe UI", 8, QFont.Normal))
        subtitle.setStyleSheet("color: #555;")
        subtitle.setAlignment(Qt.AlignCenter)

        ### Image ###
        image_label = QLabel()
        pix = QPixmap("protein.png")
        image_label.setPixmap(pix.scaled(220, 220, Qt.KeepAspectRatio, Qt.SmoothTransformation))
        image_label.setAlignment(Qt.AlignCenter)

        ### Input ###
        self.input = QLineEdit()
        self.input.setPlaceholderText("Enter Protein Uniprot ID")
        self.input.setFixedWidth(320)
        self.input.setStyleSheet("""
            QLineEdit {
                padding: 8px;
                border-radius: 8px;
                border: 2px solid #aacdf2;
                font-size: 14px;
            }
            QLineEdit:focus {
                border: 2px solid #5aa3e8;
            }
        """)

        instruction = QLabel("Example: P04637, P38398, Q9Y6K9")
        instruction.setFont(QFont("Segoe UI", 9))
        instruction.setStyleSheet("color: #666;")
        instruction.setAlignment(Qt.AlignCenter)

        ### Continue Button ###
        continue_btn = QPushButton("Continue")
        continue_btn.setFixedSize(170, 42)
        continue_btn.setStyleSheet("""
            QPushButton {
                background-color: #5aa3e8;
                color: white;
                padding: 8px;
                border-radius: 8px;
                font-size: 15px;
            }
            QPushButton:hover {
                background-color: #3a89d6;
            }
            QPushButton:pressed {
                background-color: #2f6fb2;
            }
        """)
        continue_btn.clicked.connect(self.go_to_menu)

        instruction = QLabel("Please input protein Uniprot ID")
        instruction.setFont(QFont("Arial", 9, QFont.Light))
        instruction.setAlignment(Qt.AlignCenter)

         # ===== Card Container =====
        card = QFrame()
        card.setStyleSheet("""
            QFrame {
                background: white;
                border-radius: 20px;
            }
        """)
        card.setFixedWidth(420)

        card.setFixedWidth(420)

        card_layout = QVBoxLayout(card)
        card_layout.setSpacing(18)
        card_layout.setContentsMargins(40, 40, 40, 40)
        card_layout.addWidget(title)
        card_layout.addWidget(subtitle)
        card_layout.addSpacing(10)
        card_layout.addWidget(image_label)
        card_layout.addWidget(self.input, alignment=Qt.AlignCenter)
        card_layout.addWidget(instruction)
        card_layout.addSpacing(10)
        card_layout.addWidget(continue_btn, alignment=Qt.AlignCenter)

        """layout = QVBoxLayout()
        layout.addStretch()
        layout.addWidget(title)
        layout.addWidget(image_label)
        layout.addWidget(self.input, alignment=Qt.AlignCenter)
        layout.addWidget(instruction, alignment = Qt.AlignCenter)
        layout.addWidget(continue_btn, alignment=Qt.AlignCenter)
        layout.addStretch()"""
         # ===== Main Layout =====
        layout = QVBoxLayout(self)
        layout.addStretch()
        layout.addWidget(card, alignment=Qt.AlignCenter)
        layout.addStretch()
        self.setLayout(layout)

    #def go_to_menu(self):
     #   protein_code = self.input.text().strip()
        
      #  summary_text = g1_protein.protein_summary(protein_code)
       # if summary_text is None:
        #    self.show_error_message("Protein not found. Please enter a valid UniProt ID.")
        #if protein_code:
         #   self.parent.protein_code = protein_code
          #  summary_text = g1_protein.protein_summary(protein_code)  # make sure this returns a string
           # self.parent.page2.update_summary(summary_text)
            #self.parent.stack.setCurrentIndex(1)
    def go_to_menu(self):
        protein_code = self.input.text().strip()
        
        if not protein_code:
            self.show_error_message("Please enter a protein ID.")
            return

        try:
            summary_text = g1_protein.protein_summary(protein_code)

            if summary_text is None:
                self.show_error_message(
                    "❌ Protein not found.\nPlease check the UniProt ID and try again."
                )
                return

            self.parent.protein_code = protein_code
            summary_text = g1_protein.protein_summary(protein_code)  # string
            self.parent.page2.update_summary(summary_text)
            self.parent.stack.setCurrentIndex(1)

        except Exception as e:
            self.show_error_message(f"Unexpected error: {e}")

    def show_error_message(self, message):
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Warning)
        msg.setWindowTitle("Error")
        msg.setText(message)
        msg.exec_()


#------------------------------------------------
#helping functions placed here
#---------------------------------------------
def fig_to_pixmap(fig):
    buf = BytesIO()
    try:
        fig.tight_layout()
    except Exception:
        pass
    # Save without bbox_inches='tight' so the figure scale/margins are preserved
    fig.savefig(buf, format="png", dpi=150)
    buf.seek(0)
    pixmap = QPixmap()
    pixmap.loadFromData(buf.getvalue(), "PNG")
    return pixmap

def create_card(text):
        card = QLabel(text)
        card.setWordWrap(True)
        card.setStyleSheet("""
            QLabel {
                background-color: #ffffff;
                border-radius: 10px;
                padding: 14px;
                font-family: "Segoe UI";
                font-size: 14px;
                color: #333;
            }
        """)
        return card

def dataframe_to_table(df):
    model = QStandardItemModel()
    model.setColumnCount(len(df.columns))
    model.setHorizontalHeaderLabels(df.columns.tolist())

    for row in df.itertuples(index=False):
        items = []
        for value in row:
            item = QStandardItem(str(value))
            item.setEditable(False)
            items.append(item)
        model.appendRow(items)

    proxy = QSortFilterProxyModel()
    proxy.setSourceModel(model)
    proxy.setFilterKeyColumn(-1)  # global filtering

    table = QTableView()
    table.setModel(proxy)
    table.setSortingEnabled(True)
    table.setAlternatingRowColors(True)
    table.setSelectionBehavior(QTableView.SelectRows)
    table.setSelectionMode(QTableView.SingleSelection)
    table.horizontalHeader().setStretchLastSection(False)
    table.horizontalHeader().setDefaultAlignment(Qt.AlignLeft)

    # Resize columns based on header text
    table.resizeColumnsToContents()
    table.verticalHeader().setVisible(False)

    # Clean modern style
    table.setStyleSheet("""
        QTableView {
            font-family: "Segoe UI";
            font-size: 13px;
            gridline-color: transparent;
            background-color: #ffffff;
        }
        QHeaderView::section {
            background-color: #f3f6fa;
            padding: 6px;
            border: none;
            font-weight: 600;
        }
        QTableView::item {
            padding: 6px;
        }
        QTableView::item:selected {
            background-color: #cce5ff;
        }
    """)

    return table, proxy

def create_filter_row(df, proxy):
    filter_widget = QWidget()
    layout = QHBoxLayout(filter_widget)
    layout.setContentsMargins(0, 0, 0, 0)

    for col in range(len(df.columns)):
        edit = QLineEdit()
        edit.setPlaceholderText(df.columns[col])
        edit.setClearButtonEnabled(True)

        def make_filter(column):
            return lambda text: proxy.setFilterKeyColumn(column) or proxy.setFilterFixedString(text)

        edit.textChanged.connect(make_filter(col))
        layout.addWidget(edit)

    return filter_widget

# -----------------------------------------------------------
# PAGE 2: MENU PAGE
# -----------------------------------------------------------

class MenuPage(QWidget):
    def __init__(self, parent):
        super().__init__()
        self.parent = parent
        # use the same subtle gradient background as WelcomePage
        self.setStyleSheet("""
            QWidget {
                background: qlineargradient(
                    x1:0, y1:0, x2:1, y2:1,
                    stop:0 #d9ecff,
                    stop:1 #f4f9ff
                );
            }
        """)

        # Card container (white rounded panel)
        card = QFrame()
        card.setStyleSheet("""
            QFrame {
                background: white;
                border-radius: 16px;
                border: 1px solid #e6eef8;
            }
        """)
        card.setFixedSize(760, 550)

        # Title
        title = QLabel("ProSeqNet Navigator")
        title.setFont(QFont("Segoe UI", 14, QFont.Normal))
        title.setAlignment(Qt.AlignCenter)

        # Summary box
        self.summary_box = QTextEdit()
        self.summary_box.setReadOnly(True)
        self.summary_box.setFixedHeight(300)
        self.summary_box.setStyleSheet("""
            QTextEdit {
                background: #fbfdff;
                border: 1px solid #eef6fb;
                border-radius: 8px;
                padding: 10px;
                font-family: 'Segoe UI';
                font-size: 14px;
                color: #222;
            }
        """)

        # Buttons
        self.structure_btn = QPushButton("Alphafold Structure")
        self.network_btn = QPushButton("Protein–Protein Interactions")
        self.sequence_btn = QPushButton("Sequence Analysis")
        self.variant_btn = QPushButton("Variant Analysis")

        for btn in [self.structure_btn, self.network_btn, self.sequence_btn, self.variant_btn]:
            btn.setFixedHeight(40)
            btn.setStyleSheet("""
                QPushButton {
                    background-color: #43a2d3;
                    color: white;
                    border-radius: 10px;
                    padding: 8px 14px;
                    font-size: 16px;
                }
                QPushButton:hover { background-color: #3a89d6; }
            """)

        # Connect actions
        self.structure_btn.clicked.connect(self.open_structure_dialog)
        self.network_btn.clicked.connect(self.open_ppi_dialog)
        self.variant_btn.clicked.connect(self.open_variant_dialog)

        # Back / small action
        back_btn = QPushButton("Back to Menu")
        back_btn.setFixedSize(140, 36)
        back_btn.setStyleSheet("""
            QPushButton { background-color: #92d4f6; color: #2b6ea3; border-radius: 8px; }
            QPushButton:hover { background-color: #e6f0fb; }
        """)
        back_btn.clicked.connect(self.go_back)

        # Layout inside card
        card_layout = QVBoxLayout(card)
        card_layout.setContentsMargins(28, 24, 28, 24)
        card_layout.setSpacing(14)
        card_layout.addWidget(title)
        card_layout.addWidget(self.summary_box)

        # button rows
        row1 = QHBoxLayout()
        row1.addWidget(self.structure_btn)
        row1.addWidget(self.network_btn)
        card_layout.addLayout(row1)

        row2 = QHBoxLayout()
        row2.addWidget(self.variant_btn)
        row2.addWidget(self.sequence_btn)
        card_layout.addLayout(row2)

        # footer with back button
        footer = QHBoxLayout()
        footer.addStretch()
        footer.addWidget(back_btn)
        card_layout.addLayout(footer)

        # Main layout centers the card
        layout = QVBoxLayout(self)
        layout.addStretch()
        layout.addWidget(card, alignment=Qt.AlignCenter)
        layout.addStretch()
        self.setLayout(layout)
    #------- Functions-----------
    def update_summary(self, text):
        self.summary_box.setText(text)
    
    def open_structure_dialog(self):
        pdb_data, summary_text = g1_protein.get_alphafold_pdb(self.parent.protein_code)

        if pdb_data:
            dialog = AlphaDialog(pdb_data, summary_text)
            dialog.exec_()
        else:
            QMessageBox.warning(self, "Error", "Could not load AlphaFold structure.")

    def open_ppi_dialog(self):
        fig, explain = g1_protein.ppi_network(self.parent.protein_code)

        if fig is None:
            QMessageBox.warning(self, "Error", "Could not create PPI network.")
            return
        dialog = PPIDialog("Protein–Protein Interaction Network", explain, fig)
        dialog.exec_()

    def open_variant_dialog(self):
        fig, summary_text, explain_text = g3_variant.Variant_analysis(self.parent.protein_code)

        if fig is None:
            QMessageBox.warning(self, "Error", "Could not perform variant analysis.")
            return

        dialog = VariantDialog(fig, summary_text, explain_text, self.parent.protein_code)
        dialog.exec_()


    def go_back(self):
        self.parent.stack.setCurrentIndex(0)

    


# -----------------------------------------------------------
# DIALOGUE 1: ALPHA FOLD STRUCTURE
# -----------------------------------------------------------

def AlphaDialog(pdb_data, summary_text, title="AlphaFold Structure"):

    dialog = QDialog()
    dialog.setWindowTitle(title)
    dialog.setMinimumSize(900, 800)

    layout = QVBoxLayout()

    summary_box = QTextEdit()
    summary_box.setReadOnly(True)
    summary_box.setText(summary_text)
    summary_box.setStyleSheet("background-color: #e7f2ff; font-size: 14px; padding: 8px;")
    layout.addWidget(summary_box)

    # interactive 3D structure
    view = py3Dmol.view(width=800, height=600)
    view.addModel(pdb_data, 'pdb')
    view.setStyle({'cartoon': {'color':'spectrum'}})
    view.zoomTo()
    html = view._make_html()

    webview = QWebEngineView()
    webview.setHtml(html)
    layout.addWidget(webview)

    dialog.setLayout(layout)
    return dialog
# -----------------------------------------------------------
# DIALOG 2: PROTEIN-PROTEIN INTERACTION NETWORK
# -----------------------------------------------------------
class PPIDialog(QDialog):
    def __init__(self, title, text_html, fig, parent=None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.resize(900, 800)

        layout = QVBoxLayout()

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)

        content_widget = QWidget()
        content_layout = QVBoxLayout()

        text_box = QTextEdit()
        text_box.setReadOnly(True)
        text_box.setHtml(text_html)

        # Simple whole-image zoom controls (treat figure as single image)
        pixmap = fig_to_pixmap(fig)
        image_label = QLabel()
        image_label.setAlignment(Qt.AlignCenter)

        # initial scaled pixmap
        init_w, init_h = 800, 600
        scaled_pix = pixmap.scaled(init_w, init_h, Qt.KeepAspectRatio, Qt.SmoothTransformation)
        image_label.setPixmap(scaled_pix)
        image_label.setFixedSize(scaled_pix.size())

        # zoom controls
        zoom_in = QPushButton("+")
        zoom_out = QPushButton("-")
        zoom_in.setFixedSize(36, 28)
        zoom_out.setFixedSize(36, 28)

        zoom_container = QHBoxLayout()
        zoom_container.addStretch()
        zoom_container.addWidget(zoom_out)
        zoom_container.addWidget(zoom_in)
        zoom_container.addStretch()

        # scale factor state
        scale = {'factor': 1.0}

        def do_zoom(factor):
            scale['factor'] *= factor
            w = max(100, int(init_w * scale['factor']))
            h = max(100, int(init_h * scale['factor']))
            img = pixmap.scaled(w, h, Qt.KeepAspectRatio, Qt.SmoothTransformation)
            image_label.setPixmap(img)
            image_label.setFixedSize(img.size())

        zoom_in.clicked.connect(lambda: do_zoom(1.25))
        zoom_out.clicked.connect(lambda: do_zoom(0.8))

        content_layout.addWidget(text_box)
        content_layout.addLayout(zoom_container)
        content_layout.addWidget(image_label)
        content_widget.setLayout(content_layout)
        scroll.setWidget(content_widget)

        layout.addWidget(scroll)
        self.setLayout(layout)

# -----------------------------------------------------------
# Base Dialogue Disease-associated Variants
# -----------------------------------------------------------
def DiseaseVariantDialog(df_table, summary_text, parent=None):
    dlg = QDialog(parent)
    dlg.setWindowTitle("Disease-associated Variants")
    dlg.resize(900, 800)

    layout = QVBoxLayout(dlg)

    # Summary card
    summary_card = create_card(summary_text)
    scroll = QScrollArea()
    scroll.setWidgetResizable(True)
    scroll.setWidget(summary_card)
    # make summary area occupy ~40% of dialog height and be scrollable
    summary_height = max(200, int(dlg.height() * 0.4))
    scroll.setFixedHeight(summary_height)
    layout.addWidget(scroll)

    # Table + filters
    table, proxy = dataframe_to_table(df_table)
    filter_row = create_filter_row(df_table, proxy)

    # Ensure table is visible and expands to fill remaining space
    table.setMinimumHeight(250)
    table.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
    layout.addWidget(filter_row)
    layout.addWidget(table)


    return dlg

# -----------------------------------------------------------
# DIALOG 3: Variant Analysis
# -----------------------------------------------------------
def VariantDialog(fig, summary_text, explain_text, uniprot_id, title="Variant Analysis"):
    dialog = QDialog()
    dialog.setWindowTitle(title)
    dialog.setMinimumSize(900, 800)

    layout = QVBoxLayout()

    summary_box = QTextEdit()
    summary_box.setReadOnly(True)
    summary_box.setText(summary_text)
    summary_box.setStyleSheet("background-color: #e7f2ff; font-size: 14px; padding: 8px;")
    layout.addWidget(summary_box)

    # Simple whole-image zoom controls (treat figure as single image)
    pixmap = fig_to_pixmap(fig)
    image_label = QLabel()
    image_label.setAlignment(Qt.AlignCenter)
    init_w, init_h = 800, 600
    scaled_pix = pixmap.scaled(init_w, init_h, Qt.KeepAspectRatio, Qt.SmoothTransformation)
    image_label.setPixmap(scaled_pix)
    image_label.setFixedSize(scaled_pix.size())

    # put image_label inside a scroll area so zoomed image can be panned
    image_scroll = QScrollArea()
    image_scroll.setWidgetResizable(False)
    image_scroll.setWidget(image_label)

    zoom_in = QPushButton("+")
    zoom_out = QPushButton("-")
    zoom_in.setFixedSize(36, 28)
    zoom_out.setFixedSize(36, 28)

    zoom_container = QHBoxLayout()
    zoom_container.addStretch()
    zoom_container.addWidget(zoom_out)
    zoom_container.addWidget(zoom_in)
    zoom_container.addStretch()

    scale = {'factor': 1.0}

    def do_zoom(factor):
        scale['factor'] *= factor
        w = max(100, int(init_w * scale['factor']))
        h = max(100, int(init_h * scale['factor']))
        img = pixmap.scaled(w, h, Qt.KeepAspectRatio, Qt.SmoothTransformation)
        image_label.setPixmap(img)
        image_label.setFixedSize(img.size())

    zoom_in.clicked.connect(lambda: do_zoom(1.25))
    zoom_out.clicked.connect(lambda: do_zoom(0.8))

    layout.addLayout(zoom_container)
    layout.addWidget(image_scroll)

    explain_box = QTextEdit()
    explain_box.setReadOnly(True)
    explain_box.setText(explain_text)
    explain_box.setStyleSheet("background-color: #e7f2ff; font-size: 14px; padding: 8px;")
    layout.addWidget(explain_box)

    # Button to open disease-associated variants table
    disease_btn = QPushButton("Disease-associated Variants")
    disease_btn.setFixedHeight(36)
    disease_btn.setStyleSheet("""
        QPushButton { background-color: #2d89ef; color: white; border-radius: 8px; padding: 6px 12px; }
        QPushButton:hover { background-color: #1c6fd6; }
    """)

    def open_disease_dialog():
        try:
            df_table, summary = g3_variant.disease_associated_variants(uniprot_id)
        except Exception as e:
            QMessageBox.warning(dialog, "Error", f"Could not get disease-associated variants:\n{e}")
            return

        if df_table is None or df_table.empty:
            QMessageBox.information(dialog, "No Data", "No disease-associated variants found for this protein.")
            return

        #try:
         #   html_table = df_table.to_html(index=False, escape=False)
        #except Exception:
         #   html_table = '<pre>' + df_table.to_string(index=False) + '</pre>'

        #content_html = f"<pre>{summary}</pre>" + html_table
        dlg = DiseaseVariantDialog(df_table, summary, parent=dialog)
        dlg.exec_()

    disease_btn.clicked.connect(open_disease_dialog)
    # place disease button below image
    layout.addWidget(disease_btn, alignment=Qt.AlignCenter)

    dialog.setLayout(layout)
    return dialog
# -----------------------------------------------------------
# PAGE 3: RESULT PAGE
# -----------------------------------------------------------
class ResultPage(QWidget):
    def __init__(self, parent):
        super().__init__()
        self.parent = parent
        self.setStyleSheet("background-color: #eef6ff;")

        title = QLabel("Analysis Result")
        title.setFont(QFont("Arial", 20, QFont.Bold))
        title.setAlignment(Qt.AlignCenter)

        self.text_box = QTextEdit()
        self.text_box.setReadOnly(True)
        self.text_box.setStyleSheet("""
            QTextEdit {
                font-size: 14px;
                padding: 10px;
                border: 2px solid #aaccee;
                border-radius: 10px;
                background: white;
            }
        """)

        back_btn = QPushButton("Back to Menu")
        back_btn.setFixedWidth(150)
        back_btn.setStyleSheet("""
            QPushButton {
                background-color: #5aa3e8;
                color: white;
                padding: 8px;
                border-radius: 8px;
                font-size: 14px;
            }
            QPushButton:hover {
                background-color: #3a89d6;
            }
        """)
        back_btn.clicked.connect(self.go_back)

        layout = QVBoxLayout()
        layout.addWidget(title)
        layout.addWidget(self.text_box)
        layout.addWidget(back_btn, alignment=Qt.AlignCenter)

        self.setLayout(layout)

    def set_text(self, text):
        self.text_box.setText(text)

    def go_back(self):
        self.parent.stack.setCurrentIndex(1)




# -----------------------------------------------------------
# MAIN APP
# -----------------------------------------------------------
class MainApp(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Protein Analysis Application")
        self.setGeometry(200, 100, 700, 550)

        self.protein_code = ""

        self.stack = QStackedWidget()

        self.page1 = WelcomePage(self)
        self.page2 = MenuPage(self)
        self.result_page = ResultPage(self)

        self.stack.addWidget(self.page1)
        self.stack.addWidget(self.page2)
        self.stack.addWidget(self.result_page)

        layout = QVBoxLayout()
        layout.addWidget(self.stack)
        self.setLayout(layout)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainApp()
    window.show()
    sys.exit(app.exec_())
