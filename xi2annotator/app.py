# Copyright (C) 2025  Technische Universitaet Berlin
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 3 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
# USA

import os

from flask import Flask
from flask_cors import CORS


def create_app():
    """
    Create the flask app.

    Currently only used for annotation.
    :return: flask app
    """
    app = Flask(__name__)

    # Load flask config
    if os.environ.get('FLASK_ENV', 'production') == 'development':
        app.config.from_object('xi2annotator.config.DevelopmentConfig')
    else:
        app.config.from_object('xi2annotator.config.ProductionConfig')
        try:
            app.config.from_envvar('XI2ANNOTATOR_SETTINGS')
        except (FileNotFoundError, RuntimeError):
            ...

    # add CORS header
    CORS(app, resources={
        r"/xiAnnotator/annotate/FULL": {
            "origins": "*",
            "headers": app.config['CORS_HEADERS']
        }
    })

    from xi2annotator import bp
    app.register_blueprint(bp)

    return app


# just for easy debug purposes
if __name__ == "__main__":
    import sys
    os.environ["XI2ANNOTATOR_DEBUG"] = "1"
    sys.path.insert(0, os.path.abspath(os.path.dirname(__file__) + '/..'))
    app = create_app()
    app.run(host="0.0.0.0", port=8084, debug=True)
