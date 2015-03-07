/** @file
    @brief Implementation

    @date 2014

    @author
    Sensics, Inc.
    <http://sensics.com/osvr>

*/

// Copyright 2014 Sensics, Inc.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

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
    // On Qt4, the showFullScreen() should come after the move/resize.
    // On Qt5, it comes before.
    this->showFullScreen();
    this->move(QPoint(screenres.x(), screenres.y()));
    this->resize(screenres.width(), screenres.height());
}

MainWindow::~MainWindow()
{
    delete ui;
}
