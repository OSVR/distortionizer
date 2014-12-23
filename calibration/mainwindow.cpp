#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QDesktopWidget>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    // Make the central widget have no margins, so we fill the whole region.
    QMainWindow::centralWidget()->layout()->setContentsMargins(0, 0, 0, 0);

    // Find out where the last screen lives.
    int last_screen = QApplication::desktop()->numScreens() - 1;
    QRect screenres = QApplication::desktop()->screenGeometry(last_screen);

    // Display full resolution on the last screen.
    this->move(QPoint(screenres.x(), screenres.y()));
    this->resize(screenres.width(), screenres.height());
    this->showFullScreen();
}

MainWindow::~MainWindow()
{
    delete ui;
}
