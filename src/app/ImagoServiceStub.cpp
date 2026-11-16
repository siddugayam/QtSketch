// macOS and Linux builds: Imago is not available, so every request returns an error.
#include "ImagoService.h"

ImagoService::ImagoService(QObject *parent) : QObject(parent) {}
ImagoService::~ImagoService() {}

void ImagoService::recognizeImage(const QUrl &) {
    emit imageRecognized(QString(), 0, QStringLiteral("Image recognition is only available in the Windows build."));
}
