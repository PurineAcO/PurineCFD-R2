"""PySide6 运行台: 填工况 -> 生成配置 -> 求解 -> 后处理."""

import json
import os
import subprocess
import sys
import time

from PySide6.QtCore import QProcess, QProcessEnvironment, QSettings, Qt, QTimer
from PySide6.QtGui import QFont, QFontDatabase, QTextCursor
from PySide6.QtWidgets import (
    QApplication, QCheckBox, QComboBox, QDialog, QDialogButtonBox, QDoubleSpinBox, QFileDialog,
    QFormLayout, QGridLayout, QGroupBox, QHBoxLayout, QLabel, QLineEdit, QMainWindow,
    QMessageBox, QPlainTextEdit, QProgressBar, QPushButton, QSpinBox, QTabWidget, QVBoxLayout,
    QWidget,
)

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from gui import farfield, runner

DEFAULT_L = {'cylinder': 1.0, 'naca0012': 1.0, 'oat15a': 0.23}
HERE = os.path.dirname(os.path.abspath(__file__))
POST_CURVES = (('history', 'Cl/Cd 历程'), ('residual', '残差历程'), ('cp_dist', 'Cp 分布'),
               ('cp_rms', 'Cp 脉动 RMS'), ('probe', '探针 Cp / PSD'),
               ('sequence', '末 6 帧纹影'), ('readme', 'README.md'))
CONSOLE_CAP = 800_000
STRUCTURED_CASES = ('cylinder', 'naca0012')
FIELD_LABELS = {
    'max_steps': '最大步数',
    'cfl': '伪时间 CFL',
    'dump_interval': '场输出间隔',
    'convergence_interval': '报告间隔',
    'threads': '线程数',
    'dt': '物理时间步',
    'inner': '内迭代上限',
    'inner_tol': '内迭代容差',
    'sweeps': '松弛扫描',
    'steady_iters': '定常预热步数',
    'steady_cfl': '定常 CFL',
    'wall_interval': '壁面输出间隔',
    'seed': '初始扰动',
}
CJK_FONT_FILES = (
    '/mnt/c/Windows/Fonts/Deng.ttf',
    '/mnt/c/Windows/Fonts/msyh.ttc',
    '/mnt/c/Windows/Fonts/simhei.ttf',
    '/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc',
    '/usr/share/fonts/truetype/wqy/wqy-microhei.ttc',
)

STYLE = ("""
QWidget { font-size: 17px; font-weight: 500; color: #1f2328; }
QMainWindow { background: #f5f6f8; }
QTabWidget::pane { background: #ffffff; border: 1px solid #e4e7eb;
                   border-radius: 8px; top: -1px; }
QTabBar::tab { padding: 10px 24px; margin-right: 3px; color: #59636e;
               border: 1px solid transparent; border-top-left-radius: 7px;
               border-top-right-radius: 7px; }
QTabBar::tab:selected { background: #ffffff; border-color: #e4e7eb; color: #2f81f7;
                        font-weight: 600; }
QTabBar::tab:hover:!selected { background: #eceef1; }
QGroupBox { border: 1px solid #e9ecef; border-radius: 8px; margin-top: 16px;
            padding: 16px 18px 18px 18px; font-weight: 600; }
QGroupBox::title { subcontrol-origin: margin; left: 12px; padding: 0 5px; color: #59636e; }
QLineEdit, QPlainTextEdit, QComboBox, QSpinBox, QDoubleSpinBox {
  background: #ffffff; border: 1px solid #d8dce1; border-radius: 6px; padding: 8px 12px;
  min-height: 26px;
}
QLineEdit:focus, QPlainTextEdit:focus, QComboBox:focus,
QSpinBox:focus, QDoubleSpinBox:focus { border: 1px solid #2f81f7; }
QLineEdit[readOnly="true"] { background: #f3f4f6; color: #59636e; }
QComboBox::drop-down { subcontrol-origin: padding; subcontrol-position: top right;
                       width: 30px; border: none; background: transparent; }
QComboBox::down-arrow { image: url(@ARROW@); width: 10px; height: 7px; }
QComboBox QAbstractItemView { background: #ffffff; border: 1px solid #d8dce1;
                              selection-background-color: #e8f1fd;
                              selection-color: #1f2328; outline: none; }
QSpinBox::up-button, QSpinBox::down-button,
QDoubleSpinBox::up-button, QDoubleSpinBox::down-button {
  width: 0px; border: none; background: transparent;
}
QPushButton { background: #ffffff; border: 1px solid #d8dce1; border-radius: 6px;
              padding: 9px 22px; min-height: 26px; }
QPushButton:hover { background: #f0f2f5; }
QPushButton:pressed { background: #e5e8ec; }
QPushButton:disabled { color: #a8b0b8; background: #f5f6f8; border-color: #e4e7eb; }
QPushButton#primary { background: #2f81f7; border-color: #2f81f7; color: #ffffff;
                      font-weight: 600; }
QPushButton#primary:hover { background: #2a74dd; }
QPushButton#primary:disabled { background: #bbd5fb; border-color: #bbd5fb; color: #ffffff; }
QPushButton#danger { color: #b42318; border-color: #eec9c5; }
QPushButton#danger:hover { background: #fdf4f3; }
QPushButton#tiny { padding: 4px 12px; min-height: 18px; }
QCheckBox { spacing: 8px; }
QCheckBox::indicator { width: 19px; height: 19px; border: 1px solid #c8ced5;
                       border-radius: 4px; background: #ffffff; }
QCheckBox::indicator:hover { border-color: #2f81f7; }
QCheckBox::indicator:checked { background: #ffffff; border-color: #2f81f7;
                               image: url(@CHECK@); }
QProgressBar { background: #ffffff; border: 1px solid #e4e7eb; border-radius: 6px;
               height: 18px; }
QProgressBar::chunk { background: #2f81f7; border-radius: 5px; }
QLabel#dim { color: #6b7480; }
QLabel#key { color: #6b7480; }
QLabel#value { color: #1f2328; font-weight: 600; }
QPlainTextEdit#console { background: #1b1f24; color: #d7dee6; border: 1px solid #2b3138;
                         border-radius: 8px; padding: 12px; }
QPlainTextEdit#console QScrollBar::handle:vertical { background: #3c434b; }
QScrollBar:vertical { background: transparent; width: 11px; margin: 2px; }
QScrollBar::handle:vertical { background: #c8ced5; border-radius: 5px; min-height: 30px; }
QScrollBar::handle:vertical:hover { background: #aab2ba; }
QScrollBar::add-line, QScrollBar::sub-line { height: 0; width: 0; }
"""
        .replace('@CHECK@', os.path.join(HERE, 'check.svg'))
        .replace('@ARROW@', os.path.join(HERE, 'down.svg')))


def ensure_cjk_font():
    """WSL 默认不带中文字体, 从 Windows 字体目录或系统 Noto/WQY 挂一个."""
    app = QApplication.instance()
    if app is None:
        return None
    for path in CJK_FONT_FILES:
        if not os.path.exists(path):
            continue
        fid = QFontDatabase.addApplicationFont(path)
        if fid == -1:
            continue
        families = QFontDatabase.applicationFontFamilies(fid)
        if families:
            font = app.font()
            font.setFamily(families[0])
            font.setPointSize(11)
            app.setFont(font)
            return families[0]
    print('警告: 未找到中文字体; 可 sudo apt install fonts-noto-cjk', file=sys.stderr)
    return None


class Spin(QSpinBox):
    def wheelEvent(self, event):
        event.ignore()


class DSpin(QDoubleSpinBox):
    def wheelEvent(self, event):
        event.ignore()


class Combo(QComboBox):
    def wheelEvent(self, event):
        event.ignore()


class ConfirmDialog(QDialog):
    def __init__(self, parent, rows):
        super().__init__(parent)
        self.setWindowTitle('确认运行')
        self.setMinimumWidth(560)
        grid = QGridLayout()
        grid.setHorizontalSpacing(18)
        grid.setVerticalSpacing(10)
        for i, (key, val) in enumerate(rows):
            k = QLabel(key)
            k.setObjectName('key')
            v = QLabel(str(val))
            v.setObjectName('value')
            v.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
            v.setWordWrap(True)
            grid.addWidget(k, i, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignTop)
            grid.addWidget(v, i, 1)
        grid.setColumnStretch(1, 1)
        buttons = QDialogButtonBox()
        ok = buttons.addButton('开始求解', QDialogButtonBox.ButtonRole.AcceptRole)
        ok.setObjectName('primary')
        buttons.addButton('取消', QDialogButtonBox.ButtonRole.RejectRole)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        box = QVBoxLayout(self)
        box.setContentsMargins(22, 20, 22, 18)
        box.setSpacing(18)
        box.addLayout(grid)
        box.addWidget(buttons)


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('PurineCFD-R2 运行台')
        self.resize(1280, 900)
        self.proc = None
        self.post_proc = None
        self.tail = None
        self.cfg_path = None
        self.run_dir = None
        self.run_number = None
        self.cfg_case = None
        self.run_done = False
        self.data = None
        self.t0 = None
        self.template_solver = {}
        self.fields = {}
        self.cjk = ensure_cjk_font()
        self.settings = QSettings('PurineCFD', 'runpad')
        self.setStyleSheet(STYLE)
        self.timer = QTimer(self)
        self.timer.setInterval(250)
        self.timer.timeout.connect(self.poll_log)
        self._build()
        self.load_case()

    def _mono_font(self, size):
        font = QFont()
        font.setFamilies(['DejaVu Sans Mono'] + ([self.cjk] if self.cjk else []))
        font.setStyleHint(QFont.StyleHint.Monospace)
        font.setPointSize(size)
        return font

    # ------------------------------------------------------------------ 界面
    def _build(self):
        self.tabs = QTabWidget()
        self.tabs.addTab(self._tab_case(), '① 算例')
        self.tabs.addTab(self._tab_flow(), '② 来流')
        self.tabs.addTab(self._tab_solver(), '③ 求解')
        self.tabs.addTab(self._tab_post(), '④ 后处理')
        self.tabs.addTab(self._tab_run(), '⑤ 运行')
        self.tabs.addTab(self._tab_console(), '⑥ 控制台')

        self.bar = QProgressBar()
        self.bar.setTextVisible(False)
        self.bar_label = QLabel('—')
        self.bar_label.setObjectName('dim')
        bar_row = QHBoxLayout()
        bar_row.setSpacing(12)
        bar_row.addWidget(self.bar, 1)
        bar_row.addWidget(self.bar_label)
        self.status = QLabel('就绪')
        self.status.setObjectName('dim')

        central = QWidget()
        v = QVBoxLayout(central)
        v.setContentsMargins(12, 12, 12, 12)
        v.setSpacing(10)
        v.addWidget(self.tabs, 1)
        v.addLayout(bar_row)
        v.addWidget(self.status)
        self.setCentralWidget(central)

    def _tab_case(self):
        w = QWidget()
        v = QVBoxLayout(w)
        v.setContentsMargins(16, 16, 16, 16)
        v.setSpacing(16)
        g = QGroupBox('算例')
        f = QFormLayout(g)
        f.setVerticalSpacing(12)
        self.case_combo = Combo()
        self.case_combo.addItems(runner.CASES)
        self.case_combo.setMaximumWidth(380)
        self.case_combo.currentTextChanged.connect(self.load_case)
        f.addRow('算例', self.case_combo)
        self.template_combo = Combo()
        self.template_combo.setMaximumWidth(380)
        self.template_combo.currentTextChanged.connect(self.load_template)
        f.addRow('配置模板', self.template_combo)
        v.addWidget(g)

        g2 = QGroupBox('文件')
        f2 = QFormLayout(g2)
        f2.setVerticalSpacing(12)
        self.mesh_edit = QLineEdit()
        f2.addRow('网格', self._browse_row(self.mesh_edit, '选择网格文件'))
        self.struct_edit = QLineEdit()
        f2.addRow('结构化邻接表', self._browse_row(self.struct_edit, '选择结构化邻接表'))
        self.out_edit = QLineEdit()
        self.out_edit.setReadOnly(True)
        self.out_edit.setPlaceholderText('开跑时自动创建')
        open_btn = QPushButton('打开')
        open_btn.clicked.connect(self.open_out)
        row = QHBoxLayout()
        row.addWidget(self.out_edit, 1)
        row.addWidget(open_btn)
        f2.addRow('输出目录', row)
        v.addWidget(g2)
        v.addStretch(1)
        return w

    def _tab_flow(self):
        w = QWidget()
        v = QVBoxLayout(w)
        v.setContentsMargins(16, 16, 16, 16)
        g = QGroupBox('来流条件')
        f = QFormLayout(g)
        f.setVerticalSpacing(12)

        def spin(lo, hi, dec, step, val):
            s = DSpin()
            s.setRange(lo, hi)
            s.setDecimals(dec)
            s.setSingleStep(step)
            s.setValue(val)
            s.setMaximumWidth(280)
            s.valueChanged.connect(self.refresh_far)
            return s

        self.ma = spin(1e-3, 0.999, 4, 0.01, 0.2)
        self.re = spin(1.0, 1e12, 0, 1000.0, 1000.0)
        self.temp = spin(50.0, 2000.0, 2, 5.0, 300.0)
        self.lref = spin(1e-4, 100.0, 4, 0.01, 1.0)
        self.alpha = spin(-89.0, 89.0, 2, 0.5, 0.0)
        f.addRow('马赫数 Ma', self.ma)
        f.addRow('雷诺数 Re', self.re)
        f.addRow('静温 T [K]', self.temp)
        f.addRow('参考长度 L [m]', self.lref)
        f.addRow('攻角 α [deg]', self.alpha)
        v.addWidget(g)

        g2 = QGroupBox('派生量')
        v2 = QVBoxLayout(g2)
        self.derived = QLabel()
        self.derived.setObjectName('value')
        self.derived.setFont(self._mono_font(12))
        self.derived.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        v2.addWidget(self.derived)
        v.addWidget(g2)
        v.addStretch(1)
        return w

    def _tab_solver(self):
        w = QWidget()
        v = QVBoxLayout(w)
        v.setContentsMargins(16, 16, 16, 16)
        v.setSpacing(16)
        self.grp_solver = QGroupBox('求解')
        self.form_solver = QFormLayout(self.grp_solver)
        self.form_solver.setVerticalSpacing(12)
        v.addWidget(self.grp_solver)
        self.grp_urans = QGroupBox('双时间步')
        self.form_urans = QFormLayout(self.grp_urans)
        self.form_urans.setVerticalSpacing(12)
        v.addWidget(self.grp_urans)
        v.addStretch(1)
        return w

    def _tab_post(self):
        w = QWidget()
        v = QVBoxLayout(w)
        v.setContentsMargins(16, 16, 16, 16)
        v.setSpacing(16)

        row = QHBoxLayout()
        row.setSpacing(16)
        g = QGroupBox('曲线 / 统计')
        gv = QVBoxLayout(g)
        self.pk = {}
        for key, text in POST_CURVES:
            c = QCheckBox(text)
            c.setChecked(True)
            self.pk[key] = c
            gv.addWidget(c)
        row.addWidget(g, 1)

        g2 = QGroupBox('云图面板')
        g2v = QVBoxLayout(g2)
        head = QHBoxLayout()
        for text, val in (('全选', True), ('全不选', False)):
            b = QPushButton(text)
            b.setObjectName('tiny')
            b.clicked.connect(lambda _=False, v=val: self.set_panels(v))
            head.addWidget(b)
        head.addStretch(1)
        g2v.addLayout(head)
        grid = QGridLayout()
        self.panels = {}
        for i, name in enumerate(runner.PANELS):
            c = QCheckBox(name)
            c.setChecked(True)
            self.panels[name] = c
            grid.addWidget(c, i % 5, i // 5)
        g2v.addLayout(grid)
        g2v.addStretch(1)
        row.addWidget(g2, 1)
        v.addLayout(row)

        g3 = QGroupBox('选项')
        f3 = QFormLayout(g3)
        f3.setVerticalSpacing(12)
        self.probe_pos = DSpin()
        self.probe_pos.setRange(0.0, 1.0)
        self.probe_pos.setDecimals(3)
        self.probe_pos.setSingleStep(0.05)
        self.probe_pos.setValue(0.45)
        self.probe_pos.setMaximumWidth(280)
        f3.addRow('探针位置 x/c', self.probe_pos)
        self.u_combo = Combo()
        self.u_combo.addItems(['按算例默认', '对称', '只正值'])
        self.u_combo.setMaximumWidth(280)
        f3.addRow('u 色标', self.u_combo)
        self.phase_combo = Combo()
        self.phase_combo.addItems(['auto', 'U', 'S'])
        self.phase_combo.setMaximumWidth(280)
        f3.addRow('历史段', self.phase_combo)
        v.addWidget(g3)
        v.addStretch(1)
        return w

    def _tab_run(self):
        w = QWidget()
        v = QVBoxLayout(w)
        v.setContentsMargins(16, 16, 16, 16)
        v.setSpacing(16)
        g = QGroupBox('求解器')
        f = QFormLayout(g)
        f.setVerticalSpacing(12)
        self.bin_edit = QLineEdit(self.settings.value('binary', runner.DEFAULT_BIN))
        f.addRow('可执行文件', self._browse_row(self.bin_edit, '选择可执行文件'))
        v.addWidget(g)

        g2 = QGroupBox('后处理')
        g2v = QVBoxLayout(g2)
        self.auto_post = QCheckBox('求解结束后自动执行')
        self.auto_post.setChecked(True)
        g2v.addWidget(self.auto_post)
        v.addWidget(g2)

        row = QHBoxLayout()
        row.setSpacing(12)
        self.btn_gen = QPushButton('生成配置 JSON')
        self.btn_gen.clicked.connect(self.on_generate)
        row.addWidget(self.btn_gen)
        row.addStretch(1)
        v.addLayout(row)

        row = QHBoxLayout()
        row.setSpacing(12)
        self.btn_run = QPushButton('开始运行')
        self.btn_run.setObjectName('primary')
        self.btn_run.clicked.connect(self.start_run)
        self.btn_post = QPushButton('仅跑后处理')
        self.btn_post.clicked.connect(lambda: self.run_post())
        self.btn_stop = QPushButton('停止')
        self.btn_stop.setObjectName('danger')
        self.btn_stop.clicked.connect(self.stop_run)
        self.btn_stop.setEnabled(False)
        for b in (self.btn_run, self.btn_post, self.btn_stop):
            row.addWidget(b)
        row.addStretch(1)
        v.addLayout(row)
        v.addStretch(1)
        return w

    def _tab_console(self):
        self.console = QPlainTextEdit()
        self.console.setObjectName('console')
        self.console.setReadOnly(True)
        self.console.setMaximumBlockCount(20000)
        self.console.setFont(self._mono_font(12))
        w = QWidget()
        v = QVBoxLayout(w)
        v.setContentsMargins(12, 12, 12, 12)
        v.addWidget(self.console)
        return w

    def _browse_row(self, edit, title):
        w = QWidget()
        h = QHBoxLayout(w)
        h.setContentsMargins(0, 0, 0, 0)
        h.setSpacing(8)
        h.addWidget(edit, 1)
        b = QPushButton('选择')
        b.clicked.connect(lambda: self._pick(edit, title))
        h.addWidget(b)
        return w

    def _pick(self, edit, title):
        path, _ = QFileDialog.getOpenFileName(self, title, os.path.dirname(edit.text()) or '.')
        if path:
            edit.setText(path)

    # ------------------------------------------------------------------ 载入
    def load_case(self):
        case = self.case_combo.currentText()
        self.template_combo.blockSignals(True)
        self.template_combo.clear()
        self.template_combo.addItems(runner.list_templates(case))
        self.template_combo.blockSignals(False)
        self.lref.setValue(DEFAULT_L.get(case, 1.0))
        if self.template_combo.count():
            self.load_template(self.template_combo.currentText())

    def load_template(self, name):
        if not name:
            return
        case = self.case_combo.currentText()
        tpl = runner.load_template(case, name)
        base = runner.case_dir(case)
        self.template_solver = tpl.get('solver', {})
        self._build_solver_fields(self.template_solver)
        far = tpl.get('farfield', {})
        mesh = tpl.get('io', {}).get('mesh', '')
        self.mesh_edit.setText(os.path.normpath(os.path.join(base, mesh)) if mesh else '')
        struct = tpl.get('io', {}).get('structured', '')
        self.struct_edit.setText(os.path.normpath(os.path.join(base, struct)) if struct else '')
        L = self.lref.value() or DEFAULT_L.get(case, 1.0)
        for widget, val in ((self.ma, far.get('Ma', 0.2)), (self.temp, far.get('T', 300.0)),
                            (self.alpha, far.get('alpha', 0.0))):
            widget.blockSignals(True)
            widget.setValue(val)
            widget.blockSignals(False)
        if 'p' in far:
            d = farfield.from_p(self.ma.value(), far['p'], self.temp.value(), L)
            self.re.blockSignals(True)
            self.re.setValue(round(d['Re']))
            self.re.blockSignals(False)
        self.refresh_far()

    def _build_solver_fields(self, solver):
        for form in (self.form_solver, self.form_urans):
            while form.rowCount():
                form.removeRow(0)
        self.fields = {}
        urans = solver.get('urans', {}) if isinstance(solver.get('urans'), dict) else {}
        for key, value in solver.items():
            if key == 'urans':
                continue
            self.fields[key] = self._add_solver_row(self.form_solver, key, value)
        for key, value in urans.items():
            self.fields['urans.' + key] = self._add_solver_row(self.form_urans, key, value)
        self.grp_urans.setVisible(bool(urans))

    def _add_solver_row(self, form, key, value):
        widget = self._make_field(value)
        widget.setMaximumWidth(280)
        label = FIELD_LABELS.get(key)
        if label:
            widget.setToolTip(key)
            form.addRow(label, widget)
        else:
            form.addRow(key, widget)
        return widget

    def _make_field(self, value):
        if isinstance(value, bool):
            w = QCheckBox()
            w.setChecked(value)
            return w
        if isinstance(value, int):
            w = Spin()
            w.setRange(0, 1_000_000_000)
            w.setValue(value)
            return w
        w = DSpin()
        tail = f'{value:.12f}'.rstrip('0').partition('.')[2]
        w.setDecimals(min(max(len(tail), 2), 12))
        w.setRange(0.0, 1e12)
        w.setValue(value)
        return w

    def refresh_far(self):
        d = farfield.from_re(self.ma.value(), self.re.value(), self.temp.value(),
                             self.lref.value())
        self.derived.setText(
            f"μ {d['mu']:.4e} Pa·s      a {d['a']:.2f} m/s      U∞ {d['U']:.2f} m/s\n"
            f"ρ∞ {d['rho']:.4e} kg/m³      p∞ {d['p']:.4f} Pa      q∞ {d['q']:.4f} Pa")

    def set_panels(self, on):
        for c in self.panels.values():
            c.setChecked(on)

    # ------------------------------------------------------------------ 取值
    def farfield_dict(self):
        d = farfield.from_re(self.ma.value(), self.re.value(), self.temp.value(),
                             self.lref.value())
        return dict(Ma=d['Ma'], T=d['T'], p=d['p'], alpha=self.alpha.value())

    def solver_dict(self):
        solver = json.loads(json.dumps(self.template_solver))
        for key, widget in self.fields.items():
            if isinstance(widget, QCheckBox):
                value = widget.isChecked()
            else:
                value = widget.value()
            if key.startswith('urans.'):
                solver.setdefault('urans', {})[key.split('.', 1)[1]] = value
            else:
                solver[key] = value
        return solver

    def threads(self):
        w = self.fields.get('threads')
        return int(w.value()) if isinstance(w, QSpinBox) else 32

    def post_options(self):
        return dict(
            panels=[k for k, c in self.panels.items() if c.isChecked()],
            readme=self.pk['readme'].isChecked(),
            probe_pos=self.probe_pos.value(),
            u=('auto', 'sym', 'pos')[self.u_combo.currentIndex()],
            phase=self.phase_combo.currentText(),
            wall=None,
            **{k: self.pk[k].isChecked() for k in
               ('history', 'residual', 'cp_dist', 'cp_rms', 'probe', 'sequence')},
        )

    def generate_config(self):
        """只根据当前界面信息写出配置文件, 不启动任何进程."""
        case = self.case_combo.currentText()
        mesh = self.mesh_edit.text().strip()
        if not os.path.isfile(mesh):
            QMessageBox.warning(self, '找不到网格文件', mesh)
            return None
        tpl = runner.load_template(case, self.template_combo.currentText())
        far = self.farfield_dict()
        if self.run_dir is None or self.cfg_case != case or self.run_done:
            self.run_number, self.run_dir = runner.make_run_dir(case, far['alpha'])
            self.cfg_case = case
            self.run_done = False
        self.cfg_path = runner.write_config(case, tpl, self.run_number, self.run_dir, mesh,
                                            self.struct_edit.text().strip(),
                                            self.solver_dict(), far)
        self.out_edit.setText(self.run_dir)
        return self.cfg_path

    def on_generate(self):
        cfg = self.generate_config()
        if cfg is None:
            return
        with open(cfg, encoding='utf-8') as f:
            self.append(f'$ 生成配置  {cfg}\n{f.read()}\n')
        if self.structured_missing():
            self.append('[提示] 该算例是结构化网格但未指定结构化邻接表: '
                        '求解不受影响, 后处理会退回非结构后端\n')
        self.tabs.setCurrentIndex(5)
        self.status.setText(f'已生成配置  {os.path.basename(cfg)}')

    def structured_missing(self):
        return (self.case_combo.currentText() in STRUCTURED_CASES
                and not self.struct_edit.text().strip())

    # ------------------------------------------------------------------ 运行
    def start_run(self):
        if self.proc is not None:
            return
        binary = self.bin_edit.text().strip()
        mesh = self.mesh_edit.text().strip()
        if not os.path.isfile(binary):
            QMessageBox.warning(self, '找不到可执行文件', binary)
            return
        if not os.path.isfile(mesh):
            QMessageBox.warning(self, '找不到网格文件', mesh)
            return
        case = self.case_combo.currentText()
        tpl = runner.load_template(case, self.template_combo.currentText())
        solver = self.solver_dict()
        far = self.farfield_dict()
        opt = self.post_options()
        chosen = [text for key, text in POST_CURVES if opt[key]]
        chosen += opt['panels']
        d = farfield.from_re(far['Ma'], self.re.value(), far['T'], self.lref.value())
        rows = [
            ('算例', f'{case}  ({self.template_combo.currentText()})'),
            ('网格', os.path.basename(mesh)),
            ('Ma / Re', f"{far['Ma']:g} / {self.re.value():.3g}"),
            ('T / L / α', f"{far['T']:g} K / {self.lref.value():g} m / {far['alpha']:g}°"),
            ('p∞ / U∞', f"{far['p']:.4f} Pa / {d['U']:.2f} m/s"),
            ('步数 / dt', f"{solver.get('max_steps', '-')} / "
                          f"{solver.get('urans', {}).get('dt', '-')}"),
            ('可执行文件', os.path.basename(binary)),
            ('线程', self.threads()),
            ('后处理', '、'.join(chosen) if chosen else '无'),
        ]
        if self.structured_missing():
            rows.append(('注意', '未指定结构化邻接表, 后处理将退回非结构后端'))
        if ConfirmDialog(self, rows).exec() != QDialog.DialogCode.Accepted:
            return

        cfg = self.generate_config()
        if cfg is None:
            return
        run_dir = self.run_dir
        self.bar.setMaximum(int(solver.get('max_steps', 0)))
        self.bar.setValue(0)
        self.bar_label.setText('—')
        self.console.clear()
        self.settings.setValue('binary', binary)

        env = QProcessEnvironment.systemEnvironment()
        env.insert('OMP_NUM_THREADS', str(self.threads()))
        self.proc = QProcess(self)
        self.proc.setWorkingDirectory(runner.ROOT)
        self.proc.setProcessEnvironment(env)
        self.proc.setProcessChannelMode(QProcess.ProcessChannelMode.SeparateChannels)
        self.proc.readyReadStandardError.connect(self.read_stderr)
        self.proc.finished.connect(self.solver_done)
        self.tail = runner.LogTail(os.path.join(run_dir, 'run.log'))
        self.t0 = time.monotonic()
        self.timer.start()
        self.btn_run.setEnabled(False)
        self.btn_stop.setEnabled(True)
        self.status.setText(f'求解中  {os.path.basename(run_dir)}')
        self.append(f'$ cd {runner.ROOT}\n$ OMP_NUM_THREADS={self.threads()} '
                    f'{binary} {cfg}\n\n')
        self.tabs.setCurrentIndex(5)
        self.proc.start(binary, [cfg])

    def stop_run(self):
        for p in (self.proc, self.post_proc):
            if p is not None:
                p.kill()

    def solver_done(self, code, _status):
        self.poll_log()
        self.timer.stop()
        wall = time.monotonic() - self.t0 if self.t0 else 0.0
        self.append(f'\n[求解器退出] exit code = {code}, 壁钟 {wall:.1f} s\n')
        self.proc = None
        self.run_done = True
        self.data = (self.cfg_path, self.run_dir)
        self.btn_run.setEnabled(True)
        self.btn_stop.setEnabled(bool(self.post_proc))
        self.status.setText(f'求解结束 (exit {code})')
        if self.auto_post.isChecked():
            self.run_post(wall if wall > 0.0 else None)

    def post_target(self, target=None):
        """选一个真正含有求解结果的目录: 优先刚跑完的, 其次当前生成的."""
        pairs = [(self.cfg_path, target)] if target else [self.data,
                                                          (self.cfg_path, self.run_dir)]
        for pair in pairs:
            if not pair:
                continue
            cfg, run_dir = pair
            if cfg and run_dir and os.path.isfile(os.path.join(run_dir, 'history.csv')):
                return cfg, run_dir
        return None, None

    def run_post(self, wall=None, target=None):
        if self.post_proc is not None:
            return
        cfg, run_dir = self.post_target(target)
        if run_dir is None:
            tried = [p[1] for p in (self.data, (None, self.run_dir)) if p and p[1]]
            QMessageBox.warning(self, '没有可后处理的求解结果',
                                '以下目录里都没有 history.csv:\n' + '\n'.join(tried))
            return
        if self.run_dir and run_dir != self.run_dir:
            self.append(f'\n[后处理] 目标改用 {os.path.relpath(run_dir, runner.ROOT)}'
                        f'（{os.path.basename(self.run_dir)} 里没有求解结果）\n')
        opt = self.post_options()
        if wall is not None:
            opt['wall'] = wall
        cmd = runner.build_post_cmd(cfg, self.case_combo.currentText(), run_dir, opt)
        self.post_proc = QProcess(self)
        self.post_proc.setWorkingDirectory(runner.ROOT)
        self.post_proc.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        self.post_proc.readyReadStandardOutput.connect(self.read_postout)
        self.post_proc.finished.connect(self.post_done)
        self.btn_post.setEnabled(False)
        self.status.setText('后处理中')
        self.append('\n$ ' + ' '.join(cmd) + '\n')
        self.tabs.setCurrentIndex(5)
        self.post_proc.start(cmd[0], cmd[1:])

    def read_postout(self):
        data = bytes(self.post_proc.readAllStandardOutput()).decode('utf-8', 'replace')
        self.append(data)

    def post_done(self, code, _status):
        self.append(f'\n[后处理退出] exit code = {code}\n')
        self.post_proc = None
        self.btn_post.setEnabled(True)
        self.btn_stop.setEnabled(self.proc is not None)
        self.status.setText(f'完成 (post exit {code})')

    def poll_log(self):
        if self.tail is None:
            return
        text = self.tail.poll()
        if not text:
            return
        self.append(text)
        for line in text.splitlines():
            parts = line.split()
            if len(parts) >= 2 and parts[0] == 'step':
                try:
                    value = int(parts[1])
                    self.bar.setValue(value)
                    self.bar_label.setText(f'step {value} / {self.bar.maximum()}')
                except ValueError:
                    pass
            if line.strip():
                self.status.setText(line.strip()[:150])

    def read_stderr(self):
        data = bytes(self.proc.readAllStandardError()).decode('utf-8', 'replace')
        if data.strip():
            self.append(data)

    def append(self, text):
        self.console.moveCursor(QTextCursor.MoveOperation.End)
        self.console.insertPlainText(text)
        self.console.moveCursor(QTextCursor.MoveOperation.End)
        doc = self.console.document()
        if doc.characterCount() > CONSOLE_CAP:
            cur = QTextCursor(doc)
            cur.setPosition(0)
            cur.setPosition(CONSOLE_CAP // 2, QTextCursor.MoveOperation.KeepAnchor)
            cur.removeSelectedText()

    def open_out(self):
        path = self.run_dir or runner.run_root()
        if not os.path.isdir(path):
            return
        if os.path.exists('/mnt/c/Windows/explorer.exe'):
            try:
                win = subprocess.run(['wslpath', '-w', path], capture_output=True, text=True,
                                     check=True).stdout.strip()
                subprocess.Popen(['/mnt/c/Windows/explorer.exe', win])
                return
            except Exception:
                pass
        subprocess.Popen(['xdg-open', path])

    def closeEvent(self, event):
        running = [name for name, p in (('求解器', self.proc), ('后处理', self.post_proc))
                   if p is not None]
        if running:
            ans = QMessageBox.question(
                self, '确认退出', '／'.join(running) + ' 正在运行, 关闭窗口会终止它。确定退出？',
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No)
            if ans != QMessageBox.StandardButton.Yes:
                event.ignore()
                return
        event.accept()


def main():
    app = QApplication(sys.argv)
    win = MainWindow()
    win.show()
    return app.exec()


if __name__ == '__main__':
    sys.exit(main())
