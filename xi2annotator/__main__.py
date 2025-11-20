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

"""
Main entry point for xi2annotator web service.
"""
import argparse
from waitress import serve


def main():
    """Run the annotation web service."""
    parser = argparse.ArgumentParser(description='XiAnnotator web service')
    parser.add_argument('--host', default='0.0.0.0', help='Host to bind to (default: 0.0.0.0)')
    parser.add_argument('--port', type=int, default=8084, help='Port to bind to (default: 8084)')
    parser.add_argument('--debug', action='store_true', help='Run in debug mode (Flask dev server)')
    args = parser.parse_args()

    from xi2annotator.app import create_app
    app = create_app()

    if args.debug:
        print(f"Starting debug server on {args.host}:{args.port}")
        app.run(host=args.host, port=args.port, debug=True)
    else:
        print(f"Starting production server on {args.host}:{args.port}")
        serve(app, host=args.host, port=args.port)


if __name__ == "__main__":
    main()
